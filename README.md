qme-ng: Quiver Mutation Explorer
================================

Welcome to qme-ng !

QME is a fast quiver manipulation program intended to:
- find the cardinality of the mutation class of a given quiver
- find the maximum length of a green suite starting from a given quiver

In order to do so, it uses various optimizations that are specific to the
problem at hand, and which result from a careful study by Grégoire Dupont and
Matthieu Pérotin.


Usage
------------------------------------------------------------------------------
qme-ng can work in two modes:
  1. Quiver mutation class cardinality
  2. Max green suite length
These two modes have specific options but both need to be provided a quiver to
start with. This can be done with the following options:

 --file: give a file as input. The file format can be:

* an incidence matrix in the classical Maple format
* a Keller app file

 --type x --size y: there some predifined types supported by qme-ng. These two
                    parameters both define the quiver type and its size.

1. Quiver Mutation class cardinality options

 --dump-class prefix: at the end of the exploration, each member of the mutation class will be dumped in a file with a name starting with "prefix"

2. Green suite length options

 --green: enter this mode

 --no-iso: do not use isomorphism discrimination. Much longer, but use virtually no memory

 --max\_depth: limits the tree exploration to the value passed as argument

 --p: declare the maximum length to be infinity if an edge multiplicity
      goes above the provided value


History
------------------------------------------------------------------------------
This project started in 2005 when G. Dupont and M. Pérotin where both Ph.D
students, in Mathematics and Computer Science respectively. G. Dupont needed a
quick tool to explore the quiver space while M. Pérotin needed a CPU intensive
application to serve as benchmark for the OS schedulers he was developping.

In 2007, Joris Calabrese, a computer engineer student made an internship
working on the program and a first parallel version of the mutation class
cardinality was obtained.

Now, in 2012, G. Dupont is in post doctoral studies and M. Pérotin is working
as a SW Architect for Bull. qme-ng has evolved as a green suite length computer
and is still in active developpment in order to solve memory usage issues, and
find new compute or memory friendly algorithms.

Contact
------------------------------------------------------------------------------
Matthieu Pérotin matthieu.perotin(a)bull.net
