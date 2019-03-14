const phylotree = require('phylotree');
const tape = require('tape');
const phylogenetic_sampling = require('./phylogenetic_sampling');

tape('return proper distances', function(test) {
  const newick = "(((A:.1,B:.2)N1:.3,C:.4)N2:.8,(D:.5,E:.6)N3:.7)N4;",
    tree = new phylotree.phylotree(newick),
    distances = phylogenetic_sampling.calculate_distances_from_node(tree, 'A');
  test.assert(Math.abs(distances.A - 0) < 1e-8);
  test.assert(Math.abs(distances.B - .3) < 1e-8);
  test.assert(Math.abs(distances.C - .8) < 1e-8);
  test.assert(Math.abs(distances.D - 2.4) < 1e-8);
  test.assert(Math.abs(distances.E - 2.5) < 1e-8);
  test.end();
});
