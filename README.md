# Electrical Circuit Generator
AMATH 352 Final Project (Spring 2023). Calculates voltages at nodes and current along connections in a given electrical circuit.

## Requirements
- Utilize schemdraw in an anaconda environment and draw a schematic diagram of the circuit
- Read and process text files in a python environment using String manipulation and typecasting
- Construct the required linear system for the circuit
- Take into account the voltage supplied at the source
- Solve for the voltages at the nodes and the currents along each edge
- Determine the overall current that flows through the wires from source to sink

## Linear Algebra Logic
A node-arc incidence matrix **(A)** is constructed, making the starting node's index -1 and the ending node's index +1.
A resistance matrix **(R)** is also constructed, expressing the resistances along each edge. \
Using these matrices, a vector expressing the voltages at each node **(v)** and a vector expressing the currents along each edge **(c)** can be calculated.

The voltage drop along each edge is equal to the current times the resistance, so 
**Av = -Rc**

The sum of the currents flowing into each node must be zero 
**A^T * c = 0**

This linear system can be represented as \
**| A . R ||v|** \
**| 0 A^T ||c| = 0**

Since the voltages of the source and sink node are known, the columns and rows containing them in **A** and **A^T** respectively can be removed. 
However, they must be transferred to the solution to the linear system, meaning the indices of the nodes connected to the source node in the vector must be changed to the source voltage.

Once the linear system is properly set it, it is easily solved with singular value decomposition to solve for the voltages and currents.

## Given Format of Text Files
m (number of nodes) \
node_id_1 x y \
node_id_2 x y \
... \
node_id_m x y \
n (number of edges) \
node_id node_id resistance_1 \
... \
node_id node_id resistance_n \
source_node_id source_voltage \
sink_node_id 
