#include "Node.h"
#include "Element.h"
#include <conio.h>
#include <iostream>
#include <Eigen/Sparse>

static const int PAS=200;
static const int PAS2=200;

int main()
{   setlocale(LC_ALL,"Rus");
 NodesList nodes;
 ElementsList elements;

  double minX = 0.0;
  double minY = 0.0;
  double mass[PAS];
  double mass2[PAS2];
  double maxX = 0;
  double maxY = 0;
  
  int nbXElems = 0;
  int nbYElems = 0;

  int size = (nbXElems + 1) * (nbYElems + 1);


 std::cout<<"Задаем область: " << std::endl;
 std::cout<<"Координата начало области ось   x : " << std::endl;
 std::cin >> minX;
 std::cout<<"\nКоордината конец области ось  x : " << std::endl;
 std::cin >> maxX;
 std::cout<<"\nКоордината начала области ось y : " << std::endl;
 std::cin >> minY;
 std::cout<<"\nКоордината конец области ось  y : " << std::endl;
 std::cin >> maxY;
 std::cout<<"\nКоличество элементов по оси x : " << std::endl;
 std::cin >> nbXElems;
 std::cout<<"\nКоличество элементов по оси y  : " << std::endl;
 std::cin >> nbYElems;
 std::cout<<"\nРазмер элементов по осям: " << std::endl;
 std::cout << (nbXElems + 1) * (nbYElems + 1); 
 
 

 
 
  for (int indY = 0; indY <= nbYElems; ++indY)
  {
     
    for (int indX = 0; indX <= nbXElems; ++indX)
	{
       nodes.push_back(
          Node(minX + indX * (maxX - minX) / nbXElems,   
               minY + indY * (maxY - minY) / nbYElems)); //hx,hy 
	    mass[indX]= minX + indX * (maxX - minX) / nbXElems;                             
       	mass2[indY]= minY + indY * (maxY - minY) / nbYElems;

	}


  }
  std::cout<< "\nhx:" << std::endl;    
  for(int i=0;i<=nbXElems;i++)
  std::cout << mass[i] << std::endl; 

  std::cout<< "\nhy:" << std::endl;    
  for(int j=0;j<=nbYElems;j++)
  std::cout << mass2[j] << std::endl; 
  
 
  
   for (int indY = 0; indY < nbYElems; ++indY)
   {
    for (int indX = 0; indX < nbXElems; ++indX)
    {
      Element elem(indY * (nbXElems + 1) + indX,
                   indY * (nbXElems + 1) + indX + 1,
                   (indY + 1) * (nbXElems + 1) + indX + 1,
                   (indY + 1) * (nbXElems + 1) + indX);
      elem.setNodes(nodes);
	 
      if (indX == 0)
      { // левая граница
		  elem.setBoundary(indY * (nbXElems + 1) + indX, (indX) * (nbXElems + 1)); // берем узлы так как граничное условие u|y=0 = x (x-узлы indX * (nbYElems + 1) + indY)
        elem.setBoundary((indY + 1) * (nbYElems + 1) + indX, 0.); //граничное условие u|y=1 = 0
		

      }
      if (indX == nbXElems - 1)
      { // правая граница
		  elem.setBoundary(indY * (nbXElems + 1) + indX, 0.);
		  elem.setBoundary((indY + 1) * (nbXElems + 1) + indX, 0.);

      }
      if (indY == 0)
      { // нижняя граница
		elem.setBoundary(indX * (nbYElems + 1) + indY, 0.);         // граничное условие u|x=0 = 0
        elem.setBoundary((indX + 1) * (nbYElems + 1) + indY, 1-indY * (nbXElems + 1) + indX); // граничное условие u|x=1 = 1-y (y-узлы indY * (nbXElems + 1) + indX)
		

	  }
      if (indY == nbYElems - 1)
      { // верхняя граница
		elem.setBoundary(indX * (nbYElems + 1) + indY, 0.);
        elem.setBoundary((indX + 1) * (nbYElems + 1) + indY, 0.);
		

      }

      elements.push_back(elem);
	 
	   }
   }





  std::vector<Eigen::Triplet<double> > leftSide;
  std::vector<Eigen::Triplet<double> > rightSide;
  leftSide.reserve(size);
  //Формирование матрицы жесткости
  std::cout << "матрица А:" << std::endl;
  for (ElementsList::iterator it = elements.begin();
       it != elements.end(); ++it)
  {
    it->calculate();
    MatrixNxN stiffness = it->stiffnessMatrix();
    int ids[ELEM_SIZE];
    it->nodesIndices(ids);
	  for (int i = 0; i < ELEM_SIZE; ++i)
    {
	    for (int j = 0; j < ELEM_SIZE; ++j)
      {
		Eigen::Triplet<double> left(ids[i], ids[j], stiffness(i, j));
        leftSide.push_back(left);
		std::cout << " " << stiffness(i,j);
						
      }
    }
	
  } 
   //Формирование матрицы правых частей
  std::cout << "\nматрица F:" << std::endl;
  for (ElementsList::iterator it = elements.begin();
       it != elements.end(); ++it)
  {
    it->calculate();
	VectorN stress = it->stressVector();
    int ids[ELEM_SIZE];
    it->nodesIndices(ids);
	  for (int i = 0; i < ELEM_SIZE; ++i)
    {
	   
	
      Eigen::Triplet<double> right(ids[i], 0, stress(i));
	  rightSide.push_back(right);
	  std::cout<<stress(i);
	  
	
  } 

  }
  
   
  //Формирование матрицы u - решение системы алгебрарических уравнений
  std::cout << std::endl;
  std::cout << "Вектор решений задачи u:" << std::endl;
  for (ElementsList::iterator it = elements.begin();
       it != elements.end(); ++it)
  {
    it->calculate();
    MatrixNxN stiffness = it->stiffnessMatrix();
	VectorN stress = it->stressVector();
	VectorN x;
	int ids[ELEM_SIZE];
    it->nodesIndices(ids);

	  for (int i = 0; i < ELEM_SIZE; ++i)
    { 
	  Eigen::Triplet<double> right(ids[i], 0, stress(i));
	  rightSide.push_back(right);
	     for (int j = 0; j < ELEM_SIZE; ++j)
	{
		
		Eigen::Triplet<double> left(ids[i], ids[j], stiffness(i, j));
        leftSide.push_back(left);
		x=stiffness.fullPivLu().solve(stress); 
		}

	  
    }
	std::cout << x << std::endl;
  } 
  



	getch();
}