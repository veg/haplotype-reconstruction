const phylotree = require('phylotree');


function calculate_distances_from_node (tree, node_name) {
  function recurse(node) {
    node.visited = true;
    const branch_length = +node.data.attribute;
    if(node.data.name == node_name) node.data.distance = 1e-12;
    if(node.parent && !node.parent.data.distance) {
      node.parent.data.distance = node.data.distance + branch_length;
    }
    if(node.parent && !node.data.distance) {
      node.data.distance = node.parent.data.distance + branch_length;
    }
    if(node.parent) recurse(node.parent)
    if(node.children) node.children.filter(node=>!node.visited).forEach(recurse);    
  }
  recurse(tree.get_node_by_name(node_name));
  const distances = {};
  tree.get_tips().forEach(node=>distances[node.data.name]=node.data.distance);
  return distances;
}

exports.calculate_distances_from_node = calculate_distances_from_node;

