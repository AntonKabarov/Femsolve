#pragma once
#include <vector>

struct Node 
{
  double x, y;

  Node(double theX, double theY)
    : x(theX), y(theY)
  {}
};

typedef std::vector<Node> NodesList;