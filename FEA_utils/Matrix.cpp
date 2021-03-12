#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Matrix.h"

Matrix::Matrix() { // Constructor for Any Matrix

}

Matrix::Matrix(unsigned int rows, unsigned int cols) { // Constructor for Any Matrix
    mrows = rows;
    mcols = cols;
    mate.resize(rows);
    for (unsigned int i = 0; i < mate.size(); i++)
    {
        mate[i].resize(cols, 0);
    }
}

int Matrix::readFile(std::string filename) // Constructor - reads matrix file
{
    // Check to see if input file exists
    std::ifstream myfile(filename);
    mcols = 0, mrows = 0;
    if (!myfile.is_open()) return 1;
    
    // Open the input file
    std::string dim;
    (getline(myfile, dim));
    std::stringstream dims(dim);
    std::vector<int> xy{ std::istream_iterator<int>(dims), std::istream_iterator<int>() };
    mcols = xy[0];

    if (xy.size() == 1)
    {
        mrows = mcols;
    }
    else if (xy.size() == 2)
    {
        mrows = xy[1];
    }
    else
    {
        return 2;
    }

    int rows = mrows, cols = mcols;
    double z;

    mate.resize(rows);
    for (int i = 0; i < rows; i++)
    {
        mate[i].resize(cols);
    }

    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            myfile >> z;
            mate[i][j] = z;
        }
    }
    myfile.close();    //close the file object. 
    return 0;
}

Matrix::Matrix(const Matrix& Mat2)
{
    mrows = Mat2.getrows();
    mcols = Mat2.getcols();
    mate = Mat2.mate;
}

Matrix::~Matrix()
{

}

Matrix Matrix::operator+(Matrix& mat2) // Matrix addition
{
    Matrix result(mcols, mrows);
    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            result(i, j) = mate[i][j] + mat2(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator-(Matrix& mat2) // Matrix subtraction
{
    Matrix result(mcols, mrows);
    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            result(i, j) = mate[i][j] - mat2(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator*(Matrix& mat2) // Matrix multiplication
{
    Matrix result(mcols, mrows);
    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            for (unsigned int k = 0; k < mrows; k++)
            {
                result(i, j) += mate[i][k] * mat2(k, j);
            }
        }
    }
    return result;
}

void Matrix::printRes() const // Print the result
{ 
    std::cout << "The result is" << std::endl;
    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            std::cout << mate[i][j] << " ";
        }
        std::cout << "" << std::endl;
    }
}

void Matrix::outRes(std::string filename) const // Write the result
{ 
    std::ofstream myfile(filename);
    myfile << "# vtk DataFile Version 2.0\n" << std::endl;

    for (unsigned int i = 0; i < mrows; i++)
    {
        for (unsigned int j = 0; j < mcols; j++)
        {
            myfile << mate[i][j] << " ";
        }
        
    }
    myfile.close();
}

double& Matrix::operator()(const unsigned int& row, const unsigned int& col) // Returns matrix value at given location
{
    return mate[row][col];
}

unsigned int Matrix::getrows() const // Return row
{
    return mrows;
}

unsigned int Matrix::getcols() const // Return column
{
    return mcols;
}

int Matrix::readSign(char sign, Matrix mat1, Matrix mat2) // Check sign and dimensions
{
    unsigned int r1 = mat1.getrows(), r2 = mat1.getrows(), c1 = mat1.getcols(), c2 = mat2.getcols();
    if (((sign == '+' || sign == '-') && ((r1 != r2) || (c1 != c2))) || ((sign == '*') && (r1 != c2)))
        {
            return 1;
        }
    return 0;
}
  