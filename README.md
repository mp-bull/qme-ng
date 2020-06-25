qme-ng: Quiver Mutation Explorer
================================

Welcome to qme-ng !

QME is a fast quiver manipulation program intended to:
- find the cardinality of the mutation class of a given quiver
- find the length of the green sequences starting from a given quiver

In order to do so, it uses various optimizations that are specific to the
problem at hand, and which result from a careful study by Grégoire Dupont and
Matthieu Pérotin.


Features
------------------------------------------------------------------------------

qme-ng can work in two modes:
  1. Quiver mutation class cardinality
  2. Max green sequences length

qme-ng:
 - uses isomorphisms discrimination in order to speed its exploration;
 - uses a fast isomorphism algorithm ([nauty](http://cs.anu.edu.au/~bdm/nauty/));
 - uses [arbitrary precision arithmetic libraries](http://gmplib.org/), and is thus not limited to your CPU registers size;
 - can read [Bernhard Keller's java application](http://www.math.jussieu.fr/~keller/quivermutation/) files as input;
 - produces exploitable outputs;
 - prints its mutations sequences in a format compatible with [Bernhard Keller's java application](http://www.math.jussieu.fr/~keller/quivermutation/);
 - contains tens of small optimizations to cut as many branches in the exploration tree as fast as possible;
 - is free software (BSD License, see the LICENSE file for more details);

Contact
------------------------------------------------------------------------------
Webpage: http://mp-bull.github.com/qme-ng/

Matthieu Pérotin matthieu.perotin(a)bull.net
