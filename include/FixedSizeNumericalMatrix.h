#ifndef MATHTOOLKIT_LINEARALGEBRA_FIXEDSIZENUMERICALMATRIX_H
#define MATHTOOLKIT_LINEARALGEBRA_FIXEDSIZENUMERICALMATRIX_H

#include <string>
#include <cmath>
#include <complex>

#include <array>
#include <initializer_list>

#include <vector>

#include "OutOfRangeException.h"

#include "FixedSizeNumericalVector.h"


namespace MathToolkit
{
	namespace LinearAlgebra
	{
		/*
		Column Major Layout for _columnSize-by-_rowSize matrix
		*/

		template <typename T, unsigned short column_size, unsigned short row_size>
		class NumericalMatrix final
		{
		public:
			using size_type = unsigned short;
			using value_type = T;
			using iterator = typename std::array<value_type, (column_size * row_size)>::iterator;
			using const_iterator = typename std::array<value_type, (column_size * row_size)>::const_iterator;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, operators, and destructor

		public:
			NumericalMatrix() noexcept;
			explicit NumericalMatrix(const value_type) noexcept;

			explicit NumericalMatrix(const std::vector<NumericalVector<T, column_size>>&);
			explicit NumericalMatrix(const NumericalVector<value_type, (column_size * row_size)>&) noexcept;
			explicit NumericalMatrix(NumericalVector<value_type, (column_size * row_size)>&&) noexcept;
			NumericalMatrix(const double scalingMultiplier, const NumericalVector<value_type, column_size>&, const NumericalVector<value_type, row_size>&) noexcept;

			NumericalMatrix& operator=(const value_type scalar) noexcept;

			NumericalMatrix(const NumericalMatrix&) = default;
			NumericalMatrix(NumericalMatrix&&)  noexcept = default;
			NumericalMatrix& operator=(const NumericalMatrix&) = default;
			NumericalMatrix& operator=(NumericalMatrix&&)  noexcept = default;

			virtual ~NumericalMatrix() = default;

			// Constructors, operators, and destructor
	// **********************************************************************************************************************************************************************************************************************************************************************************************
	// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Access methods

			value_type operator()(const size_type rowIndex, const size_type columnIndex) const noexcept;
			value_type& operator()(const size_type rowIndex, const size_type columnIndex) noexcept;
			value_type at(const size_type rowIndex, const size_type columnIndex) const;
			value_type& at(const size_type rowIndex, const size_type columnIndex);

			iterator begin() noexcept;
			const_iterator begin() const noexcept;
			iterator end() noexcept;
			const_iterator end() const noexcept;

			size_type columnSize() const noexcept;
			size_type rowSize() const noexcept;

		// Access methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Unary operators

			NumericalMatrix& operator+=(const NumericalMatrix&) noexcept;
			NumericalMatrix& operator-=(const NumericalMatrix&) noexcept;
			NumericalMatrix& operator*=(const NumericalMatrix&) noexcept;
			NumericalMatrix& operator*=(const value_type scalingMultiplier) noexcept;
			NumericalMatrix& operator/=(const value_type scalingDivider) noexcept;

		// Unary operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Mathematical methods

			NumericalMatrix& add(const value_type scalingMultiplier, const NumericalMatrix&) noexcept;
			NumericalMatrix& add(const value_type scalingMultiplier, const NumericalMatrix&, const value_type selfScalingMultiplier) noexcept;
			NumericalMatrix& transform(const NumericalMatrix&) noexcept;
			NumericalMatrix& transform(const value_type scalingMultiplier, const NumericalMatrix&) noexcept;
			NumericalMatrix& multiply(const value_type scalingMultiplier, const NumericalMatrix&) noexcept;

			NumericalMatrix getTransposedMatrix() const noexcept;

		// Mathematical methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Utility

			NumericalVector<T, column_size> columnVector(const size_type columnIndex) const;
			NumericalVector<T, row_size> rowVector(const size_type rowIndex) const;
			NumericalVector<T, (column_size * row_size)> toNumericalVector() noexcept;

			bool isSquareMatrix() const noexcept;
			bool isSame(const NumericalMatrix&, const double valuePrecision = std::numeric_limits<double>::epsilon()) const noexcept;

			value_type getTrace() const;
			double getHilbertSchmidtNorm() const noexcept;
			double getHilbertSchmidtNormSquare() const noexcept;

			void swap(NumericalMatrix&) noexcept;

		// Utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			std::array<value_type, (column_size * row_size)> _valueArray;
		};



		template <typename T, unsigned short rank>
		inline NumericalMatrix<T, rank, rank> getIdentityMatrix() noexcept
		{
			NumericalMatrix<T, rank, rank> identityMatrix;
			{
				for (unsigned short index = 0; index < rank; ++index)
					identityMatrix(index, index) = T{ 1.0 };
			}

			return identityMatrix;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> operator+(NumericalMatrix<T, column_size, row_size> nma, const NumericalMatrix<T, column_size, row_size>& nmb) noexcept
		{
			nma += nmb;
			return nma;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> operator-(NumericalMatrix<T, column_size, row_size> nma, const NumericalMatrix<T, column_size, row_size>& nmb) noexcept
		{
			nma -= nmb;
			return nma;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> operator*(const double scalingMultiplier, NumericalMatrix<T, column_size, row_size> nm) noexcept
		{
			nm *= scalingMultiplier;
			return nm;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalVector<T, column_size> operator*(const NumericalMatrix<T, column_size, row_size>& nm, NumericalVector<T, row_size> nv) noexcept
		{
			NumericalVector<T, column_size> resultVector(0.0);
			{
				auto matrixIter = nm.begin();
				{
					for (const auto& val : nv)
					{
						auto newIter = resultVector.begin();

						for (unsigned short rowIndex = 0; rowIndex < column_size; ++rowIndex)
						{
							(*newIter) += (*matrixIter) * val;
							++newIter;
							++matrixIter;
						}
					}
				}
			}

			return resultVector;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> operator*(NumericalMatrix<T, column_size, row_size> nma, const NumericalMatrix<T, column_size, row_size>& nmb) noexcept
		{
			nma *= nmb;
			return nma;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> operator/(NumericalMatrix<T, column_size, row_size> nm, const double scalingDivider) noexcept
		{
			nm /= scalingDivider;
			return nm;
		}

		template <typename T, unsigned short column_size, unsigned short row_size>
		inline NumericalMatrix<T, column_size, row_size> multiply(const double scalingMultiplier, NumericalMatrix<T, column_size, row_size> nma, const NumericalMatrix<T, column_size, row_size>& nmb) noexcept
		{
			nma.multiply(scalingMultiplier, nmb);
			return nma;
		}
	}
}

// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors and operators

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix() noexcept
	: _valueArray{}
{
	_valueArray.fill(value_type{ 0.0 });
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix(const value_type value) noexcept
	: _valueArray{}
{
	_valueArray.fill(value);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix(const std::vector<NumericalVector<T, column_size>>& numericalVectors)
	: _valueArray{}
{
	if(numericalVectors.size() != row_size)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "The sizes of \"numericalVectors\" is not equal to row size." };

	else
	{
		auto thisIter = std::begin(_valueArray);

		for (const auto& numericalVector : numericalVectors)
		{
			for (const auto value : numericalVector)
			{
				(*thisIter) = value;
				++thisIter;
			}
		}
	}
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix(const NumericalVector<T, (column_size * row_size)>& numericalVector) noexcept
	: _valueArray{}
{
	_valueArray = numericalVector._valueArray;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix(NumericalVector<T, (column_size * row_size)>&& numericalVector) noexcept
	: _valueArray{}
{
	_valueArray = std::move(numericalVector._valueArray);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::NumericalMatrix(const double scalingMultiplier, const NumericalVector<value_type, column_size>& nvba, const NumericalVector<value_type, row_size>& nvbb) noexcept
	: _valueArray{}
{
	auto thisIter = std::begin(_valueArray);

	for (auto rightVectorIter = nvbb.begin(); rightVectorIter != nvbb.end(); ++rightVectorIter)
	{
		for (auto leftVectorIter = nvba.begin(); leftVectorIter != nvba.end(); ++leftVectorIter)
		{
			(*thisIter) = (*leftVectorIter) * (*rightVectorIter);
			++thisIter;
		}
	}
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator=(const value_type scalar) noexcept
{
	_valueArray.fill(scalar);
	return *this;
}

// Constructors and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Access methods

template <typename T, unsigned short column_size, unsigned short row_size>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator()(const size_type rowIndex, const size_type columnIndex) const noexcept
{
	return _valueArray[(columnIndex * column_size) + rowIndex];
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline T& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator()(const size_type rowIndex, const size_type columnIndex) noexcept
{
	return _valueArray[(columnIndex * column_size) + rowIndex];
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::at(const size_type rowIndex, const size_type columnIndex) const
{
	if (rowIndex < column_size && columnIndex < row_size)
		return _valueArray[(columnIndex * column_size) + rowIndex];
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "at", "\"index\" is out of range." };
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline T& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::at(const size_type rowIndex, const size_type columnIndex)
{
	if (rowIndex < column_size && columnIndex < row_size)
		return _valueArray[(columnIndex * column_size) + rowIndex];
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "at", "\"index\" is out of range." };
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::begin() noexcept
{
	return std::begin(_valueArray);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::const_iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::begin() const noexcept
{
	return std::begin(_valueArray);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::end() noexcept
{
	return std::end(_valueArray);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::const_iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::end() const noexcept
{
	return std::end(_valueArray);
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::size_type MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::columnSize() const noexcept
{
	return column_size;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::size_type MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::rowSize() const noexcept
{
	return row_size;
}

// Access methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Unary operators

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator+=(const NumericalMatrix& numericalMatrix) noexcept
{
	for (size_type index = 0; index < (column_size * row_size); ++index)
		_valueArray[index] += numericalMatrix._valueArray[index];

	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator-=(const NumericalMatrix& numericalMatrix) noexcept
{
	for (size_type index = 0; index < (column_size * row_size); ++index)
		_valueArray[index] -= numericalMatrix._valueArray[index];

	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator*=(const NumericalMatrix& nm) noexcept
{
	NumericalMatrix resultMatrix;
	{
		for (size_type rowIndex = 0; rowIndex < resultMatrix.columnSize(); ++rowIndex)
		{
			auto nmIter = nm.begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;


			for (size_type columnIndex = 0; (1 + columnIndex) < resultMatrix.rowSize(); ++columnIndex)
			{
				auto thisIter = begin();
				thisIter += rowIndex;

				for (size_type index = 0; (1 + index) < row_size; ++index)
				{
					(*resultIter) += (*thisIter) * (*nmIter);
					thisIter += column_size;
					++nmIter;
				}

				(*resultIter) += (*thisIter) * (*nmIter);
				resultIter += resultMatrix.columnSize();
			}


			auto thisIter = begin();
			thisIter += rowIndex;

			for (size_type index = 0; (1 + index) < row_size; ++index)
			{
				(*resultIter) += (*thisIter) * (*nmIter);
				thisIter += column_size;
				++nmIter;
			}

			(*resultIter) += (*thisIter) * (*nmIter);
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator*=(const value_type scalingMultiplier) noexcept
{
	for (size_type index = 0; index < (column_size * row_size); ++index)
		_valueArray[index] *= scalingMultiplier;

	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::operator/=(const value_type scalingDivider) noexcept
{
	const value_type scalingMultiplier = (1.0 / scalingDivider);

	for (size_type index = 0; index < (column_size * row_size); ++index)
		_valueArray[index] *= scalingMultiplier;

	return *this;
}

// Unary operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Mathematical methods

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::add(const value_type scalingMultiplier, const NumericalMatrix& numericalMatrix) noexcept
{
	for (size_type index = 0; index < (column_size * row_size); ++index)
		_valueArray[index] += (numericalMatrix._valueArray[index] * scalingMultiplier);

	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::add(const value_type scalingMultiplier, const NumericalMatrix& numericalMatrix, const value_type selfScalingMultiplier) noexcept
{
	for (size_type index = 0; index < (column_size * row_size); ++index)
	{
		_valueArray[index] *= selfScalingMultiplier;
		_valueArray[index] += (numericalMatrix._valueArray[index] * scalingMultiplier);
	}

	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::transform(const NumericalMatrix& nm) noexcept
{
	NumericalMatrix resultMatrix;
	{
		for (std::size_t rowIndex = 0; rowIndex < resultMatrix.columnSize(); ++rowIndex)
		{
			auto thisIter = begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;

			for (std::size_t columnIndex = 0; (1 + columnIndex) < resultMatrix.rowSize(); ++columnIndex)
			{
				auto nmIter = nm.begin();
				nmIter += rowIndex;

				for (std::size_t index = 0; (1 + index) < nm.rowSize(); ++index)
				{
					(*resultIter) += (*nmIter) * (*thisIter);
					nmIter += nm.columnSize();
					++thisIter;
				}

				(*resultIter) += (*nmIter) * (*thisIter);
				resultIter += resultMatrix.columnSize();
			}


			auto nmIter = nm.begin();
			nmIter += rowIndex;

			for (std::size_t index = 0; (1 + index) < nm.rowSize(); ++index)
			{
				(*resultIter) += (*nmIter) * (*thisIter);
				nmIter += nm.columnSize();
				++thisIter;
			}

			(*resultIter) += (*nmIter) * (*thisIter);
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::transform(const value_type scalingMultiplier, const NumericalMatrix& nm) noexcept
{
	NumericalMatrix resultMatrix;
	{
		for (std::size_t rowIndex = 0; rowIndex < resultMatrix.columnSize(); ++rowIndex)
		{
			auto thisIter = begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;

			for (std::size_t columnIndex = 0; (1 + columnIndex) < resultMatrix.rowSize(); ++columnIndex)
			{
				auto nmIter = nm.begin();
				nmIter += rowIndex;

				for (std::size_t index = 0; (1 + index) < nm.rowSize(); ++index)
				{
					(*resultIter) += (*nmIter) * (*thisIter);
					nmIter += nm.columnSize();
					++thisIter;
				}

				(*resultIter) += (*nmIter) * (*thisIter);
				(*resultIter) *= scalingMultiplier;
				resultIter += resultMatrix.columnSize();
			}


			auto nmIter = nm.begin();
			nmIter += rowIndex;

			for (std::size_t index = 0; (1 + index) < nm.rowSize(); ++index)
			{
				(*resultIter) += (*nmIter) * (*thisIter);
				nmIter += nm.columnSize();
				++thisIter;
			}

			(*resultIter) += (*nmIter) * (*thisIter);
			(*resultIter) *= scalingMultiplier;
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>& MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::multiply(const value_type scalingMultiplier, const NumericalMatrix& nm) noexcept
{
	NumericalMatrix resultMatrix;
	{
		for (size_type rowIndex = 0; rowIndex < resultMatrix.columnSize(); ++rowIndex)
		{
			auto nmIter = nm.begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;


			for (size_type columnIndex = 0; (1 + columnIndex) < resultMatrix.rowSize(); ++columnIndex)
			{
				auto thisIter = begin();
				thisIter += rowIndex;

				for (size_type index = 0; index < row_size; ++index)
				{
					(*resultIter) += (*thisIter) * (*nmIter);
					thisIter += column_size;
					++nmIter;
				}

				(*resultIter) += (*thisIter) * (*nmIter);
				(*resultIter) *= scalingMultiplier;
				resultIter += resultMatrix.columnSize();
			}


			auto thisIter = begin();
			thisIter += rowIndex;

			for (size_type index = 0; index < row_size; ++index)
			{
				(*resultIter) += (*thisIter) * (*nmIter);
				thisIter += column_size;
				++nmIter;
			}

			(*resultIter) += (*thisIter) * (*nmIter);
			(*resultIter) *= scalingMultiplier;
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size> MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::getTransposedMatrix() const noexcept
{
	NumericalMatrix<T, row_size, column_size> transposedMatrix;
	{
		auto thisIter = std::begin(_valueArray);

		for (std::size_t columnIndex = 0; columnIndex < row_size; ++columnIndex)
		{
			auto transposedIter = transposedMatrix.begin();
			transposedIter += columnIndex;

			for (std::size_t rowIndex = 0; (1 + rowIndex) < column_size; ++rowIndex)
			{
				(*transposedIter) = (*thisIter);
				transposedIter += row_size;
				++thisIter;
			}

			(*transposedIter) = (*thisIter);
			++thisIter;
		}
	}

	return transposedMatrix;
}

// Mathematical methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Utility

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalVector<T, column_size> MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::columnVector(const size_type columnIndex) const
{
	if (columnIndex < row_size)
	{
		NumericalVector<T, column_size> columnVector;
		{
			auto iter = _valueArray.begin();
			iter += (column_size * columnIndex);

			for (size_type columnIndex = 0; columnIndex < column_size; ++columnIndex)
			{
				columnVector[columnIndex] = (*iter);
				++iter;
			}
		}

		return columnVector;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "operator[]", "\"columnIndex\" is out of range." };
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalVector<T, row_size> MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::rowVector(const size_type rowIndex) const
{
	if (rowIndex < column_size)
	{
		NumericalVector<T, row_size> rowVector;
		{
			auto iter = _valueArray.begin();
			iter += rowIndex;

			for (size_type rowIndex = 0; (1 + rowIndex) < row_size; ++rowIndex)
			{
				rowVector[rowIndex] = (*iter);
				iter += column_size;
			}

			rowVector[rowIndex] = (*iter);
		}

		return rowVector;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "rowVector", "\"rowIndex\" is out of range." };
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline MathToolkit::LinearAlgebra::NumericalVector<T, (column_size* row_size)> MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::toNumericalVector() noexcept
{
	NumericalVector<T, (column_size * row_size)> numericalVector;
	numericalVector._valueArray = std::move(_valueArray);
	return numericalVector;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::isSquareMatrix() const noexcept
{
	if (column_size == row_size)
		return true;
	else
		return false;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::isSame(const NumericalMatrix& numericalMatrix, const double valuePrecision) const noexcept
{
	auto otherIter = numericalMatrix.begin();

	for (auto thisIter = begin(); thisIter != end(); ++thisIter)
	{
		double differenceSquare = 0.0;
		{
			double differenceAbs = std::abs((*otherIter) - (*thisIter));
			differenceSquare += (differenceAbs * differenceAbs);
		}

		if (differenceSquare < (valuePrecision * valuePrecision))
			++otherIter;
		else
			return false;
	}

	return true;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::getTrace() const
{
	if (isSquareMatrix())
	{
		value_type value = 0.0;
		{
			for (size_type rowIndex = 0; rowIndex < column_size; ++rowIndex)
				value += operator()(rowIndex, rowIndex);
		}

		return value;
	}

	else
		throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "getTrace", "This is not a square matrix." };
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline double MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::getHilbertSchmidtNorm() const noexcept
{
	return std::sqrt(getHilbertSchmidtNormSquare());
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline double MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::getHilbertSchmidtNormSquare() const noexcept
{
	double normSquare{ 0.0 };
	{
		for (const auto& value : _valueArray)
		{
			double absValue = std::abs(value);
			normSquare += (absValue * absValue);
		}
	}

	return normSquare;
}

template <typename T, unsigned short column_size, unsigned short row_size>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, column_size, row_size>::swap(NumericalMatrix& numericalMatrix) noexcept
{
	_valueArray.swap(numericalMatrix._valueArray);
}

// Utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************

#endif // !MATHTOOLKIT_LINEARALGEBRA_FIXEDSIZENUMERICALMATRIX_H
