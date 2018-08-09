#pragma once
#include "Node.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

static const int ELEM_SIZE = 4;

typedef Eigen::Matrix<double, 2, ELEM_SIZE> Matrix2xN;
typedef Eigen::Matrix<double, ELEM_SIZE, 2> MatrixNx2;
typedef Eigen::Matrix<double, ELEM_SIZE, ELEM_SIZE> MatrixNxN;
typedef Eigen::Matrix<double, ELEM_SIZE, 1> VectorN;
class Element
{
private:
  int nodesID[ELEM_SIZE];
  Matrix2xN nodes;
  MatrixNxN stiffness;
  VectorN stress;

  bool isBound[ELEM_SIZE];
  double boundValue[ELEM_SIZE];

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Element(int theNode1, int theNode2,
          int theNode3, int theNode4);

  void setNodes(const NodesList& theAllNodes);






  // ���������� ������� ���������
void nodesIndices(int theIDs[ELEM_SIZE]);

  // ���������� ������� �������� � ������� ��������
  void calculate();

  // ��������� ������� ��������
  MatrixNxN stiffnessMatrix();

  // ��������� ������ ��������
  VectorN stressVector();

  // ������� ���������� ������� � ����
  void setBoundary(int nodeIndex, double value);

private:
  // ������� �����
  VectorN formFunctions(double xi, double eta) const;

  // �������� ������� �����
  MatrixNx2 formFunctionDerivatives(double xi, double eta) const;

  // ����� ������
  Eigen::Matrix2Xd gaussPoints(int order = 2) const;

  // ���� ������
  Eigen::VectorXd gaussWeights(int order = 2) const;

  double sigma(double x, double y) const
  { return 1.; }
  double f(double x, double y) const
  { return 0.; }
};

typedef std::vector<Element, Eigen::aligned_allocator<Element>  > ElementsList;