#include "Element.h"

Element::Element(int theNode1, int theNode2,
                 int theNode3, int theNode4)
{
	//theNode1 - x1 theNode2 - y1 theNode3 - x2 theNode4 - y2
  nodesID[0] = theNode1;
  nodesID[1] = theNode2;
  nodesID[2] = theNode3;
  nodesID[3] = theNode4;

  isBound[0] = false;
  isBound[1] = false;
  isBound[2] = false;
  isBound[3] = false;

  boundValue[0] = 0.; // U|x=0 = 0
  boundValue[1] = 1-theNode2; // U|x=1 = 1-y (y - theNode2)
  boundValue[2] = theNode1; // u|y=0 = 0
  boundValue[3] = 0.; // u|y=1 = 0
}
void Element::nodesIndices(int theIDs[ELEM_SIZE])
{
	for(int i=0;i<ELEM_SIZE;i++)
	{  theIDs[i]=i;
	
	}	

}

void Element::setNodes(const NodesList& theAllNodes)
{
  for (int i = 0; i < ELEM_SIZE; ++i)
  {
    nodes.coeffRef(0, i) = theAllNodes[nodesID[i]].x;
    nodes.coeffRef(1, i) = theAllNodes[nodesID[i]].y;
  }
}

void Element::setBoundary(int nodeIndex, double value)
{
  for (int i = 0; i < ELEM_SIZE; ++i)
    if (nodesID[i] == nodeIndex)
    {
      isBound[i] = true;
      boundValue[i] = value;
    }
}

MatrixNxN Element::stiffnessMatrix()
{
  return stiffness;
}

VectorN Element::stressVector()
{
  return stress;
}

void Element::calculate()
{
  Eigen::Matrix2Xd gaussP = gaussPoints();
  Eigen::VectorXd  gaussW = gaussWeights();

  stiffness.setZero();
  stress.setZero();

  for (int p = 0; p < gaussP.rows(); ++p)
  {
    double xi = gaussP(0, p);
    double eta = gaussP(1, p);
    double w = gaussW(p);
	
    VectorN N = formFunctions(xi, eta);
    MatrixNx2 dN = formFunctionDerivatives(xi, eta);

    Eigen::Matrix2d J = nodes * dN;
    double detJ = J.determinant();
    Eigen::Matrix2d invJ = J.inverse();

    Eigen::Vector2d xy = nodes * N;
    stiffness.noalias() += sigma(xy(0), xy(1)) * dN * invJ * invJ.transpose() * dN.transpose() * detJ * w;
    stress.noalias() += f(xy(0), xy(1)) * N * detJ * w;
	
  }

  // применяем граничные условия
  for (int bnd = 0; bnd < ELEM_SIZE; ++bnd)
  {
    if (isBound[bnd])
    {
      stiffness.row(bnd).setZero();
      stiffness.coeffRef(bnd, bnd) = 1.0;
      stress.coeffRef(bnd) = boundValue[bnd];
    }
  }
}

VectorN Element::formFunctions(double xi, double eta) const
{
  VectorN func;
  func << (1.-xi)/2. * (1.-eta)/2.,
          (1.+xi)/2. * (1.-eta)/2.,
          (1.+xi)/2. * (1.+eta)/2.,
          (1.-xi)/2. * (1.+eta)/2.;
  return func;
}

MatrixNx2 Element::formFunctionDerivatives(double xi, double eta) const
{
  MatrixNx2 deriv;
  deriv << -0.25 * (1.-eta), -0.25 * (1.-xi),
            0.25 * (1.-eta), -0.25 * (1.+xi),
            0.25 * (1.+eta),  0.25 * (1.+xi),
           -0.25 * (1.+eta),  0.25 * (1.-xi);
  return deriv;
}

Eigen::Matrix2Xd Element::gaussPoints(int order) const
{
  switch (order)
  {
  case 1: {
    Eigen::Matrix2Xd gauss(2, 1);
    gauss << 0., 0.;
    return gauss;
          }
  case 2: {
    Eigen::Matrix2Xd gauss(2, 4);
    gauss << -1./sqrt(3.), -1./sqrt(3.),
              1./sqrt(3.), -1./sqrt(3.),
              1./sqrt(3.),  1./sqrt(3.),
             -1./sqrt(3.),  1./sqrt(3.);
    return gauss;
          }
  }
  return Eigen::Matrix2Xd();
}


Eigen::VectorXd Element::gaussWeights(int order) const
{
  switch (order)
  {
  case 1: {
    Eigen::VectorXd gauss(1);
    gauss << 4.;
    return gauss;
          }
  case 2: {
    Eigen::VectorXd gauss(4);
    gauss << 1., 1., 1., 1.;
    return gauss;
          }
  }
  return Eigen::VectorXd();
}
