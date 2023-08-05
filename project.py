import numpy as np
from scipy.linalg import lu_factor, lu_solve
import schemdraw
import schemdraw.elements as elm
from schemdraw import flow
schemdraw.use('svg')

file_name = "test-data/test1.txt"
data = np.array(open(file_name, "r").readlines()) # Data imported into a numpy array

num_nodes = int(data[0].replace("\n","")) # The first value is the number of nodes
coord = {} # Stores all of the coordinates in a dictionary sorted by the node_id
node_ids = [] # List so index of the node_id can be accessed in the dictionary

# Removes the spaces in the array
for i in range(1, num_nodes+1):
    coord.update({data[i].split(" ")[0]: data[i].split(" ")[1:]}) 
    node_ids.append(data[i].split(" ")[0]) # The node_id is appended to the node_ids list


# Removes the "\n" and type-casts all values to ints
for i in coord:
    for j in range(len(coord.get(i))):
        coord.get(i)[j] = float(coord.get(i)[j].replace("\n",""))
print("What the coordinate dictionary looks like currently: \n" + str(coord))
data = np.delete(data, np.arange(0, num_nodes + 1)) # Removes processed data from the array
num_edges = int(data[0]) # The first remaining value is the number of edges
data = np.delete(data, 0) # Since it has been processed, it's safe to be removed
a = np.zeros((num_edges, num_nodes)) # The node-arc incidence matrix
r = np.diagflat(np.ones(num_edges)) # The resistance matrix
source_connect = {} # Keeps track of where the source connects to update the rhs later

v0 = [data[num_edges].split(" ")[0], int(data[num_edges].split(" ")[-1])] # Source node
sink = data[num_edges+1].split(" ")[0].replace("\n", "") # Sink node
for i in range(0, num_edges):
    curr = data[i].split(" ") # Creates a list that can be more easily accessed
    start = int(node_ids.index(curr[0])) # Start node_id index of the edge
    end = int(node_ids.index(curr[1])) # End node_id index of the edge
    
    a[i, start] = -1 # Updates start node in the node-arc array
    a[i, end] = 1 # Updates end node in the node-arc array
    
    
    if (curr[0] == v0[0]): # If the source node is the start of the current edge
        source_connect.update({i : v0[1]})
    if (curr[1] == v0[0]): # If the source node is the end of the current edge
        source_connect.update({i : -v0[1]})
    
    resistance = int(curr[2].replace("\n","")) # resistance of the current edge
    r[i] = r[i] * resistance # Updates the edge in the resistance array

print("Node-arc Incidence Matrix: \n"+str(a)+"\n")
print("Resistance Matrix: \n"+str(r)+"\n")
print("Source-Connection Dictionary: "+str(source_connect))
# Checking the Source and Sink variables:
print("Source Index:",v0[0],"at coordinates", coord.get(v0[0]), "with", v0[1], "volts")
print("Sink Index:", sink, "at coordinates", coord.get(sink))

schematic = schemdraw.Drawing(canvas='svg') # Instantiates the drawing
ynew = 4
if (coord.get(v0[0])[1] < 0):
    ynew = -4
schematic += elm.SourceV().endpoints(coord.get(v0[0]), (coord.get(v0[0])[0], 
coord.get(v0[0])[1]+ynew)).label(str(v0[1])+'V') # Creates the source at the v0 node

for i in coord.keys():
    schematic += elm.Line().endpoints(coord.get(i), 
    coord.get(i)).label(i, color = "#48cae4", loc = "center") 
    # Creates a node at each coordinate
schematic += elm.Ground().at(coord.get(sink)) # Creates the sink node at the given node_id

for i in range(num_edges):
    index=[coord.get(node_ids[list(a[i]).index(-1)]),coord.get(node_ids[list(a[i]).index(1)])]
    schematic += (R1 := elm.Resistor().endpoints(index[0], index[1])) 
    # Resistance is found at the diagonal of the given indices
    schematic += elm.CurrentLabelInline(direction='in').at(R1).label('$'+
            str(int(r[i,i]))+'\Omega$', rotate = True , color = "#023e8a")
    # Resistance is found at the diagonal of the given 
schematic += flow.Start().label('Original Plot').at([1.5, 2.5])
schematic.draw()

a1 = np.delete(a, [node_ids.index(v0[0]), node_ids.index(sink)], axis = 1)
arr1 = np.concatenate([a1,r], axis = 1, dtype = float) # Top row of the array is concatenated
arr2 = np.concatenate([np.zeros((num_nodes-2, num_nodes-2)), a1.T], axis = 1, dtype = float) 
# Bottom row of the array is concatenated
arr = np.concatenate([arr1, arr2], axis = 0, dtype = float) # Top and bottom rows concatenated
rhs = np.zeros(num_nodes + num_edges - 2) # Solution is 0s
start_node = [node_ids.index(v0[0])]
# Fix the rhs using the source_connect dictionary
for i in source_connect:
    rhs[i] = source_connect.get(i)
# Solution not using a source_connect dictionary:
    # rhs[:4] = rhs[:4] - (a[:, node_ids.index(v0[0])]) * v0[1]
print("A: \n",arr)  
print("Right-Hand Side:\n", rhs)
print("Source Connect Dictionary: "+str(source_connect))
print("Updated Zero Vector: "+str(rhs))

# LU Factorization
lu, piv = lu_factor(arr)
avec = lu_solve((lu, piv), rhs)
print("-"*20, "\nLU Factorization Solution:\n", avec)

# Singular Value Decomposition
u, s, vt = np.linalg.svd(arr)
utb = u.T @ rhs / s
soln = sum(utb[i] * vt[i] for i in range(len(utb)))
print("SVD Solution:", soln)

def format(arr_in):
    # Formatting the voltages
    if (node_ids.index(v0[0]) < node_ids.index(sink)):
        voltages = arr_in[:node_ids.index(v0[0])] # Everything before the source node
        voltages = np.append(voltages, v0[1]) # Source node
        voltages = np.append(voltages, arr_in[node_ids.index(v0[0]):node_ids.index(sink)-1]) 
        # Between
        voltages = np.append(voltages, 0) # Sink node
        voltages = np.append(voltages, arr_in[node_ids.index(sink)-1:num_nodes-2]) 
        # Everything before the sink node
    else:
        voltages = arr_in[:node_ids.index(sink)] # Everything before the source node
        voltages = np.append(voltages, 0) # Source node
        voltages = np.append(voltages, soln[node_ids.index(sink):node_ids.index(v0[0])-1]) 
        # Between
        voltages = np.append(voltages, v0[1]) # Sink node
        voltages = np.append(voltages, arr_in[node_ids.index(v0[0])-1:num_nodes-2]) 
        # Everything before the sink node

    currents = arr_in[num_nodes-2:] # Currents
    return voltages, currents

v_lu, c_lu = format(avec)
v_svd, c_svd = format(soln)

# Printing the x vector to check
print("Voltages")
for i in range(num_nodes):
    print(node_ids[i], str(v_lu[i]), "volts")
print("\nCurrents")
for i in range(num_edges):
    print(node_ids[list(a[i]).index(-1)], "-->", node_ids[list(a[i]).index(1)], c_lu[i], "amps")

# Definition of the Schematic Drawing
def draw_schem(vers, schematic, voltages, currents, method):
    ynew = 4
    if (coord.get(v0[0])[1] < 0):
        ynew = -4
    schematic += elm.SourceV().endpoints(coord.get(v0[0]), (coord.get(v0[0])[0], 
    coord.get(v0[0])[1]+ynew)).label(str(v0[1])+'V') 
    # Creates the source at the v0 node
    for i in coord.keys(): # Node IDs
        if vers.lower() == "voltage":
            schematic += elm.Line().endpoints(coord.get(i), coord.get(i)).label(i+
            " - "+str(voltages[node_ids.index(i)])[:3]+"V", color = "#48cae4", loc = "center")
            schematic += flow.Start(w=5, color = "#023e8a").label(method+' Voltage Drawing').at([1.5, 2.5])
        else :
            schematic += elm.Line().endpoints(coord.get(i), coord.get(i)).label(i,
                                                color = "#48cae4", loc = "center") 
            schematic += flow.Start(w=5, color = "#0096c7").label(method+' Current Drawing').at([1.5, 2.5])
        # Creates a node at each coordinate
        schematic += elm.Ground().at(coord.get(sink)) 
        # Creates the sink node at the given node_id
    
    for i in range(num_edges): 
        index=[coord.get(node_ids[list(a[i]).index(-1)]),
                coord.get(node_ids[list(a[i]).index(1)])]
            # Index of edges found by which indices of the node-arc array are -1 and 1
        if vers.lower() == "current":
            schematic += (R1 := elm.Resistor().endpoints(index[0], index[1])).label('$'+
            str(int(r[i,i]))+'\Omega$', loc = "bottom", color = "#023e8a") 
            # Resistance is found at the diagonal of the given indices
            schematic += elm.CurrentLabelInline(direction='in').at(R1).label(
                str(currents[i])[:5], color = "#0096c7", rotate = True)
        else:
            schematic += (R1 := elm.Resistor().endpoints(index[0], index[1]))

    return schematic

draw_schem("voltage", schemdraw.Drawing(canvas='svg'), v_lu, c_lu, "LU").draw()
draw_schem("current", schemdraw.Drawing(canvas='svg'), v_lu, c_lu, "LU").draw()

draw_schem("voltage", schemdraw.Drawing(canvas='svg'), v_svd, c_svd, "SVD").draw()
draw_schem("current", schemdraw.Drawing(canvas='svg'), v_svd, c_svd, "SVD").draw()

print("\n", "-"*20, "\nLU Factorization Verification")
# The voltage drop along each edge is equal to the current times the resistance
print("Predicted:", str(a @ v_lu))
print("Expected:", str(-r @ c_lu))
print("Residuals:", str(a @ v_lu + r @ c_lu), "Sum (should be 0):", (a @ v_lu + r @ c_lu).sum())
# The sum of the currents flowing into each node must be zero
print("Residuals: " + str(a1.T @ c_lu), "Sum (should be 0):", (a1.T @ c_lu).sum())

print("\n", "-"*20, "\nSVD Verification")
# The voltage drop along each edge is equal to the current times the resistance
print("Predicted:", str(a @ v_svd))
print("Expected:", str(-r @ c_svd))
print("Residuals:", str(a @ v_svd + r @ c_svd), "Sum (should be 0):", (a @ v_svd + r @ c_svd).sum())
# The sum of the currents flowing into each node must be zero
print("Residuals:", str(a1.T @ c_svd), "Sum (should be 0):", (a1.T @ c_svd).sum())