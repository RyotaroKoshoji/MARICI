#ifndef MATHTOOLKIT_LINEARALGEBRA_NUMERICALMATRIX_H
#define MATHTOOLKIT_LINEARALGEBRA_NUMERICALMATRIX_H

#include <string>
#include <type_traits>
#include <initializer_list>

#include <algorithm>
#include <limits>
#include <complex>
#include <valarray>
#include <cmath>

#include <utility>
#include <vector>

#include "ArgumentOutOfRangeException.h"
#include "InvalidOperationException.h"

#include "NumericalVector.h"
#include "NumericalVectorSlice.h"

#include "FixedSizeNumericalMatrix.h"


namespace MathToolkit
{
	namespace LinearAlgebra
	{
		/*
		Column Major Layout for _columnSize-by-_rowSize matrix
		*/

		template <typename T>
		class NumericalMatrix<T, 0, 0> final
		{
		public:
			using size_type = std::size_t;
			using value_type = T;
			using iterator = value_type*;
			using const_iterator = const value_type*;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, operators, and destructor

		public:
			NumericalMatrix() noexcept;
			NumericalMatrix(const size_type columnSize, const size_type rowSize);
			NumericalMatrix(const size_type columnSize, const size_type rowSize, const value_type);

			explicit NumericalMatrix(const std::vector<NumericalVector<value_type, 0>>&);
			NumericalMatrix(const double scalingMultiplier, const NumericalVector<value_type, 0>&, const NumericalVector<value_type, 0>&);
			NumericalMatrix(const size_type columnSize, const size_type rowSize, const NumericalVector<value_type, 0>&);
			NumericalMatrix(const size_type columnSize, const size_type rowSize, NumericalVector<value_type, 0>&&);

			NumericalMatrix(const NumericalMatrix&) = default;
			NumericalMatrix(NumericalMatrix&&) noexcept = default;
			NumericalMatrix& operator=(const NumericalMatrix&) = default;
			NumericalMatrix& operator=(NumericalMatrix&&) noexcept = default;

			template <typename U>
			NumericalMatrix(const NumericalMatrix<U, 0, 0>);

			template <typename U>
			NumericalMatrix& operator=(const NumericalMatrix<U, 0, 0>);

			NumericalMatrix& operator=(const value_type scalar) noexcept;


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

			NumericalVectorSlice<value_type> columnVectorSlice(const size_type columnIndex);
			NumericalVectorSlice<value_type> rowVectorSlice(const size_type rowIndex);

		// Access methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Unary operators

			NumericalMatrix& operator+=(const NumericalMatrix&);
			NumericalMatrix& operator-=(const NumericalMatrix&);
			NumericalMatrix& operator*=(const NumericalMatrix&);
			NumericalMatrix& operator*=(const value_type scalingMultiplier) noexcept;
			NumericalMatrix& operator/=(const value_type scalingDivider) noexcept;

		// Unary operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Mathematical methods

			NumericalMatrix& add(const value_type scalingMultiplier, const NumericalMatrix&);
			NumericalMatrix& add(const value_type scalingMultiplier, const NumericalMatrix&, const value_type selfScalingMultiplier);
			NumericalMatrix& transform(const NumericalMatrix&);
			NumericalMatrix& transform(const value_type scalingMultiplier, const NumericalMatrix&);
			NumericalMatrix& multiply(const value_type scalingMultiplier, const NumericalMatrix&);

			NumericalMatrix& invert();
			NumericalMatrix getInverseMatrix() const;
			NumericalMatrix getTransposedMatrix() const noexcept;

			void getSingularValueDecomposition() const;

		// Mathematical methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Utility

			NumericalVector<value_type, 0> columnVector(const size_type columnIndex) const;
			NumericalVector<value_type, 0> rowVector(const size_type rowIndex) const;
			NumericalVector<value_type, 0> toNumericalVector() noexcept;

			bool isEmpty() const noexcept;
			bool isSquareMatrix() const noexcept;
			bool isSame(const NumericalMatrix&, const double valuePrecision = std::numeric_limits<double>::epsilon()) const noexcept;

			value_type getTrace() const;
			double getHilbertSchmidtNorm() const noexcept;
			double getHilbertSchmidtNormSquare() const noexcept;

		// Utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Customize methods

			void swap(NumericalMatrix&) noexcept;
			void resize(const size_type columnSize, const size_type rowSize, const value_type value = value_type{});
			void clear() noexcept;

			void pushbackColumn(const NumericalVector<value_type, 0>&);
			void pushbackRow(const NumericalVector<value_type, 0>&);
			void insertColumn(const NumericalVector<value_type, 0>&, const size_type columnIndex);
			void insertRow(const NumericalVector<value_type, 0>&, const size_type rowIndex);
			void eraseColumn(const size_type columnIndex);
			void eraseRow(const size_type rowIndex);

		// Customize methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			std::valarray<value_type> _valueArray;

			size_type _columnSize;
			size_type _rowSize;
		};


		template <typename T>
		inline NumericalMatrix<T, 0, 0> getIdentityMatrix(const std::size_t rank) noexcept
		{
			NumericalMatrix<T, 0, 0> identityMatrix(rank, rank);
			{
				for (std::size_t index = 0; index < rank; ++index)
					identityMatrix(index, index) = T{ 1.0 };
			}

			return identityMatrix;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> operator+(NumericalMatrix<T, 0, 0> nma, const NumericalMatrix<T, 0, 0>& nmb)
		{
			nma += nmb;
			return nma;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> operator-(NumericalMatrix<T, 0, 0> nma, const NumericalMatrix<T, 0, 0>& nmb)
		{
			nma -= nmb;
			return nma;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> operator*(const double scalingMultiplier, NumericalMatrix<T, 0, 0> nm) noexcept
		{
			nm *= scalingMultiplier;
			return nm;
		}

		template <typename T>
		inline NumericalVector<T, 0> operator*(const NumericalMatrix<T, 0, 0>& nm, NumericalVector<T, 0> nv)
		{
			nv.transform(nm);
			return nv;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> operator*(NumericalMatrix<T, 0, 0> nma, const NumericalMatrix<T, 0, 0>& nmb)
		{
			nma *= nmb;
			return nma;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> operator/(NumericalMatrix<T, 0, 0> nm, const double scalingDivider) noexcept
		{
			nm /= scalingDivider;
			return nm;
		}

		template <typename T>
		inline NumericalMatrix<T, 0, 0> multiply(const double scalingMultiplier, NumericalMatrix<T, 0, 0> nma, const NumericalMatrix<T, 0, 0>& nmb)
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

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix() noexcept
	: _valueArray{}
	, _columnSize{ 0 }
	, _rowSize{ 0 }
{
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const size_type columnSize, const size_type rowSize)
	: _valueArray(value_type{ 0.0 }, (columnSize * rowSize))
	, _columnSize{ columnSize }
	, _rowSize{ rowSize }
{
	if (_columnSize == 0 || _rowSize == 0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "\"columnSize\" and/or \"rowSize\" are zero." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const size_type columnSize, const size_type rowSize, const value_type value)
	: _valueArray(value, (columnSize * rowSize))
	, _columnSize{ columnSize }
	, _rowSize{ rowSize }
{
	if (_columnSize == 0 || _rowSize == 0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "\"columnSize\" and/or \"rowSize\" are zero." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const std::vector<NumericalVector<T, 0>>& numericalVectors)
	: _valueArray(value_type{ 0.0 }, (numericalVectors.at(0).size() * numericalVectors.size()))
	, _columnSize{ numericalVectors.at(0).size() }
	, _rowSize{ numericalVectors.size() }
{
	size_type columnSize = numericalVectors[0].size();
	auto thisIter = std::begin(_valueArray);

	for (const auto& numericalVector : numericalVectors)
	{
		if (columnSize == numericalVector.size())
		{
			for (const auto value : numericalVector)
			{
				(*thisIter) = value;
				++thisIter;
			}
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "The column size of \"numericalVectors\" is not consistent." };
	}
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const size_type columnSize, const size_type rowSize, const NumericalVector<T, 0>& numericalVector)
	: _valueArray(value_type{ 0.0 }, (columnSize* rowSize))
	, _columnSize{ columnSize }
	, _rowSize{ rowSize }
{
	if (_valueArray.size() == numericalVector.size())
		_valueArray = numericalVector._valueArray;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "The size of \"numericalVector\" is not equal to (\"columnSize\" * \"rowSize\")." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const size_type columnSize, const size_type rowSize, NumericalVector<T, 0>&& numericalVector)
	: _valueArray(value_type{ 0.0 }, (columnSize * rowSize))
	, _columnSize{ columnSize }
	, _rowSize{ rowSize }
{
	if (_valueArray.size() == numericalVector.size())
		_valueArray = std::move(numericalVector._valueArray);
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "constructor", "The size of \"numericalVector\" is not equal to (\"columnSize\" * \"rowSize\")." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const double scalingMultiplier, const NumericalVector<value_type, 0>& nvba, const NumericalVector<value_type, 0>& nvbb)
	: _valueArray(value_type{ 0.0 }, (nvba.size()* nvbb.size()))
	, _columnSize{ nvba.size() }
	, _rowSize{ nvbb.size() }
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

template <typename T>
template <typename U>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::NumericalMatrix(const NumericalMatrix<U, 0, 0> nm)
	: _valueArray(value_type{ 0.0 }, (nm.columnSize()* nm.rowSize()))
	, _columnSize{ nm.columnSize() }
	, _rowSize{ nm.rowSize() }
{
	auto thisIter = std::begin(_valueArray);

	for (auto inputIter = nm.begin(); inputIter != nm.end(); ++inputIter, ++thisIter)
		*thisIter = static_cast<T>(*inputIter);
}

template <typename T>
template <typename U>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator=(const NumericalMatrix<U, 0, 0> nm)
{
	_valueArray.resize((nm.columnSize() * nm.rowSize()), value_type{});
	_columnSize = nm.columnSize();
	_rowSize = nm.rowSize();
	{
		auto thisIter = std::begin(_valueArray);

		for (auto inputIter = nm.begin(); inputIter != nm.end(); ++inputIter, ++thisIter)
			*thisIter = static_cast<T>(*inputIter);
	}

	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator=(const value_type scalar) noexcept
{
	_valueArray = scalar;
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

template <typename T>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator()(const size_type rowIndex, const size_type columnIndex) const noexcept
{
	return _valueArray[(columnIndex * _columnSize) + rowIndex];
}

template <typename T>
inline T& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator()(const size_type rowIndex, const size_type columnIndex) noexcept
{
	return _valueArray[(columnIndex * _columnSize) + rowIndex];
}

template <typename T>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::at(const size_type rowIndex, const size_type columnIndex) const
{
	if (rowIndex < _columnSize && columnIndex < _rowSize)
		return _valueArray[(columnIndex * _columnSize) + rowIndex];
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "at", "\"index\" is out of range." };
}

template <typename T>
inline T& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::at(const size_type rowIndex, const size_type columnIndex)
{
	if (rowIndex < _columnSize && columnIndex < _rowSize)
		return _valueArray[(columnIndex * _columnSize) + rowIndex];
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "at", "\"index\" is out of range." };
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::begin() noexcept
{
	return std::begin(_valueArray);
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::const_iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::begin() const noexcept
{
	return std::begin(_valueArray);
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::end() noexcept
{
	return std::end(_valueArray);
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::const_iterator MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::end() const noexcept
{
	return std::end(_valueArray);
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::size_type MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::columnSize() const noexcept
{
	return _columnSize;
}

template <typename T>
inline typename MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::size_type MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::rowSize() const noexcept
{
	return _rowSize;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalVectorSlice<T> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::columnVectorSlice(const size_type columnIndex)
{
	if (columnIndex < _rowSize)
		return NumericalVectorSlice(_valueArray, std::slice(columnIndex * _columnSize, _columnSize, 1));
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "columnVectorSlice", "\"columnIndex\" is out of range." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalVectorSlice<T> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::rowVectorSlice(const size_type rowIndex)
{
	if (rowIndex < _columnSize)
		return NumericalVectorSlice(_valueArray, std::slice(rowIndex, _rowSize, _columnSize));
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "rowVectorSlice", "\"rowIndex\" is out of range." };
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

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator+=(const NumericalMatrix& numericalMatrix)
{
	_valueArray += numericalMatrix._valueArray;
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator-=(const NumericalMatrix& numericalMatrix)
{
	_valueArray -= numericalMatrix._valueArray;
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator*=(const NumericalMatrix& nm)
{
	NumericalMatrix resultMatrix(_columnSize, nm.rowSize(), 0.0);
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

				for (size_type index = 0; (1 + index) < _rowSize; ++index)
				{
					(*resultIter) += (*thisIter) * (*nmIter);
					thisIter += _columnSize;
					++nmIter;
				}

				(*resultIter) += (*thisIter) * (*nmIter);
				resultIter += resultMatrix.columnSize();
			}


			auto thisIter = begin();
			thisIter += rowIndex;

			for (size_type index = 0; (1 + index) < _rowSize; ++index)
			{
				(*resultIter) += (*thisIter) * (*nmIter);
				thisIter += _columnSize;
				++nmIter;
			}


			(*resultIter) += (*thisIter) * (*nmIter);
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator*=(const value_type scalingMultiplier) noexcept
{
	_valueArray *= scalingMultiplier;
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::operator/=(const value_type scalingDivider) noexcept
{
	const value_type scalingMultiplier = (1.0 / scalingDivider);
	_valueArray *= scalingMultiplier;
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

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::add(const value_type scalingMultiplier, const NumericalMatrix& numericalMatrix)
{
	auto iter = std::begin(numericalMatrix._valueArray);

	for (auto& value : _valueArray)
	{
		value += (*iter) * scalingMultiplier;
		++iter;
	}

	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::add(const value_type scalingMultiplier, const NumericalMatrix& numericalMatrix, const value_type selfScalingMultiplier)
{
	auto iter = std::begin(numericalMatrix._valueArray);

	for (auto& value : _valueArray)
	{
		value *= selfScalingMultiplier;
		value += (*iter) * scalingMultiplier;
		++iter;
	}

	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::transform(const NumericalMatrix& nm)
{
	NumericalMatrix resultMatrix(nm._columnSize, _rowSize);
	{
		for (std::size_t rowIndex = 0; rowIndex < resultMatrix._columnSize; ++rowIndex)
		{
			auto thisIter = begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;

			for (std::size_t columnIndex = 0; (1 + columnIndex) < resultMatrix._rowSize; ++columnIndex)
			{
				auto nmIter = nm.begin();
				nmIter += rowIndex;

				for (std::size_t index = 0; (1 + index) < nm._rowSize; ++index)
				{
					(*resultIter) += (*nmIter) * (*thisIter);
					nmIter += nm._columnSize;
					++thisIter;
				}

				(*resultIter) += (*nmIter) * (*thisIter);
				resultIter += resultMatrix._columnSize;
			}


			auto nmIter = nm.begin();
			nmIter += rowIndex;

			for (std::size_t index = 0; (1 + index) < nm._rowSize; ++index)
			{
				(*resultIter) += (*nmIter) * (*thisIter);
				nmIter += nm._columnSize;
				++thisIter;
			}

			(*resultIter) += (*nmIter) * (*thisIter);
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::transform(const value_type scalingMultiplier, const NumericalMatrix& nm)
{
	NumericalMatrix resultMatrix(nm._columnSize, _rowSize);
	{
		for (std::size_t rowIndex = 0; rowIndex < resultMatrix._columnSize; ++rowIndex)
		{
			auto thisIter = begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;

			for (std::size_t columnIndex = 0; (1 + columnIndex) < resultMatrix._rowSize; ++columnIndex)
			{
				auto nmIter = nm.begin();
				nmIter += rowIndex;

				for (std::size_t index = 0; (1 + index) < nm._rowSize; ++index)
				{
					(*resultIter) += (*nmIter) * (*thisIter);
					nmIter += nm._columnSize;
					++thisIter;
				}

				(*resultIter) += (*nmIter) * (*thisIter);
				(*resultIter) *= scalingMultiplier;
				resultIter += resultMatrix._columnSize;
			}


			auto nmIter = nm.begin();
			nmIter += rowIndex;

			for (std::size_t index = 0; (1 + index) < nm._rowSize; ++index)
			{
				(*resultIter) += (*nmIter) * (*thisIter);
				nmIter += nm._columnSize;
				++thisIter;
			}

			(*resultIter) += (*nmIter) * (*thisIter);
			(*resultIter) *= scalingMultiplier;
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>& MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::multiply(const value_type scalingMultiplier, const NumericalMatrix& nm)
{
	NumericalMatrix resultMatrix(_columnSize, nm._rowSize, 0.0);
	{
		for (size_type rowIndex = 0; rowIndex < resultMatrix._columnSize; ++rowIndex)
		{
			auto nmIter = nm.begin();
			auto resultIter = resultMatrix.begin();
			resultIter += rowIndex;


			for (size_type columnIndex = 0; (1 + columnIndex) < resultMatrix._rowSize; ++columnIndex)
			{
				auto thisIter = begin();
				thisIter += rowIndex;

				for (size_type index = 0; (1 + index) < _rowSize; ++index)
				{
					(*resultIter) += (*thisIter) * (*nmIter);
					thisIter += _columnSize;
					++nmIter;
				}

				(*resultIter) += (*thisIter) * (*nmIter);
				(*resultIter) *= scalingMultiplier;
				resultIter += resultMatrix._columnSize;
			}


			auto thisIter = begin();
			thisIter += rowIndex;

			for (size_type index = 0; (1 + index) < _rowSize; ++index)
			{
				(*resultIter) += (*thisIter) * (*nmIter);
				thisIter += _columnSize;
				++nmIter;
			}

			(*resultIter) += (*thisIter) * (*nmIter);
			(*resultIter) *= scalingMultiplier;
		}
	}

	this->swap(resultMatrix);
	return *this;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::getInverseMatrix() const
{
	NumericalMatrix inverseMatrix{ *this };
	inverseMatrix.invert();
	return inverseMatrix;
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::getTransposedMatrix() const noexcept
{
	NumericalMatrix<T, 0, 0> transposedMatrix(_rowSize, _columnSize);
	{
		auto thisIter = std::begin(_valueArray);

		for (std::size_t columnIndex = 0; columnIndex < _rowSize; ++columnIndex)
		{
			auto transposedIter = transposedMatrix.begin();
			transposedIter += columnIndex;

			for (std::size_t rowIndex = 0; (1 + rowIndex) < _columnSize; ++rowIndex)
			{
				(*transposedIter) = (*thisIter);
				transposedIter += _rowSize;
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

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalVector<T, 0> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::columnVector(const size_type columnIndex) const
{
	if (columnIndex < _rowSize)
	{
		NumericalVector<T, 0> columnVector;
		columnVector._valueArray = _valueArray[std::slice(columnIndex * _columnSize, _columnSize, 1)];
		return columnVector;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "operator[]", "\"columnIndex\" is out of range." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalVector<T, 0> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::rowVector(const size_type rowIndex) const
{
	if (rowIndex < _columnSize)
	{
		NumericalVector<T, 0> rowVector;
		rowVector._valueArray = _valueArray[std::slice(rowIndex, _rowSize, _columnSize)];
		return rowVector;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "rowVector", "\"rowIndex\" is out of range." };
}

template <typename T>
inline MathToolkit::LinearAlgebra::NumericalVector<T, 0> MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::toNumericalVector() noexcept
{
	NumericalVector<T, 0> numericalVector;
	numericalVector._valueArray = std::move(_valueArray);
	return numericalVector;
}

template <typename T>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::isEmpty() const noexcept
{
	return !(static_cast<bool>(_valueArray.size()));
}

template <typename T>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::isSquareMatrix() const noexcept
{
	if (isEmpty())
		return false;
	else
	{
		if (_columnSize == _rowSize)
			return true;
		else
			return false;
	}
}

template <>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<double, 0, 0>::isSame(const NumericalMatrix& numericalMatrixBase, const double valuePrecision) const noexcept
{
	if ((_columnSize == numericalMatrixBase._columnSize) && (_rowSize == numericalMatrixBase._rowSize))
	{
		auto otherIter = numericalMatrixBase.begin();

		for (auto thisIter = begin(); thisIter != end(); ++thisIter)
		{
			double differenceSquare = 0.0;
			{
				value_type difference = ((*otherIter) - (*thisIter));
				differenceSquare += (difference * difference);
			}

			if (differenceSquare < (valuePrecision * valuePrecision))
				++otherIter;
			else
				return false;
		}

		return true;
	}

	else
		return false;
}

template <>
inline bool MathToolkit::LinearAlgebra::NumericalMatrix<std::complex<double>, 0, 0>::isSame(const NumericalMatrix& numericalMatrixBase, const double valuePrecision) const noexcept
{
	if ((_columnSize == numericalMatrixBase._columnSize) && (_rowSize == numericalMatrixBase._rowSize))
	{
		auto otherIter = numericalMatrixBase.begin();

		for (auto thisIter = begin(); thisIter != end(); ++thisIter)
		{
			double differenceSquare = 0.0;
			{
				value_type difference = ((*otherIter) - (*thisIter));
				differenceSquare += (difference.real() * difference.real());
				differenceSquare += (difference.imag() * difference.imag());
			}

			if (differenceSquare < (valuePrecision * valuePrecision))
				++otherIter;
			else
				return false;
		}

		return true;
	}

	else
		return false;
}


template <typename T>
inline T MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::getTrace() const
{
	if (isSquareMatrix())
	{
		value_type value = 0.0;
		{
			for (size_type rowIndex = 0; rowIndex < _columnSize; ++rowIndex)
				value += operator()(rowIndex, rowIndex);
		}

		return value;
	}

	else
		throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "getTrace", "This is not a square matrix." };
}

template <typename T>
inline double MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::getHilbertSchmidtNorm() const noexcept
{
	return std::sqrt(getHilbertSchmidtNormSquare());
}

template <>
inline double MathToolkit::LinearAlgebra::NumericalMatrix<double, 0, 0>::getHilbertSchmidtNormSquare() const noexcept
{
	double normSquare{ 0.0 };
	{
		for (const auto value : _valueArray)
			normSquare += (value * value);
	}

	return normSquare;
}

template <>
inline double MathToolkit::LinearAlgebra::NumericalMatrix<std::complex<double>, 0, 0>::getHilbertSchmidtNormSquare() const noexcept
{
	double normSquare{ 0.0 };
	{
		for (const auto value : _valueArray)
		{
			normSquare += (value.real() * value.real());
			normSquare += (value.imag() * value.imag());
		}
	}

	return normSquare;
}

// Utility
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
// Customize methods

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::swap(NumericalMatrix& nm) noexcept
{
	std::swap(_columnSize, nm._columnSize);
	std::swap(_rowSize, nm._rowSize);
	_valueArray.swap(nm._valueArray);
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::resize(const size_type columnSize, const size_type rowSize, const value_type value)
{
	_columnSize = columnSize;
	_rowSize = rowSize;
	_valueArray.resize((_columnSize * _rowSize), value);
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::clear() noexcept
{
	std::valarray<value_type> va;
	_valueArray.swap(va);
	_columnSize = 0;
	_rowSize = 0;
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::pushbackColumn(const NumericalVector<T, 0>& nv)
{
	insertColumn(nv, rowSize());
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::pushbackRow(const NumericalVector<T, 0>& nv)
{
	insertRow(nv, columnSize());
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::insertColumn(const NumericalVector<T, 0>& columnVector, const size_type columnIndex)
{
	if (columnIndex <= _rowSize)
	{
		if (isEmpty())
		{
			_valueArray = columnVector._valueArray;
			_columnSize = columnVector.size();
			_rowSize = 1;
		}

		else
		{
			std::valarray<value_type> newValueArray(_columnSize * (_rowSize + 1));
			{
				auto insertIter = std::begin(_valueArray);
				insertIter += (_columnSize * columnIndex);

				auto newInsertIter = std::begin(newValueArray);
				newInsertIter += (_columnSize * columnIndex);

				std::copy(std::begin(_valueArray), insertIter, std::begin(newValueArray));
				std::copy(columnVector.begin(), columnVector.end(), newInsertIter);
				std::copy(insertIter, std::end(_valueArray), (newInsertIter += _columnSize));
			}

			_valueArray.swap(newValueArray);
			_rowSize += 1;
		}
	}


	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "insertColumn", "\"columnIndex\" is out of range." };
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::insertRow(const NumericalVector<T, 0>& rowVector, const size_type rowIndex)
{
	if (rowIndex <= _columnSize)
	{
		if (isEmpty())
		{
			_valueArray = rowVector._valueArray;
			_columnSize = 1;
			_rowSize = rowVector.size();
		}

		else
		{
			std::valarray<value_type> newValueArray((_columnSize + 1) * _rowSize);
			{
				auto iter = std::begin(_valueArray);
				auto newIter = std::begin(newValueArray);
				auto rowVectorIter = rowVector.begin();
				size_type remainedSize = _columnSize - rowIndex;

				for (size_type columnIndex = 0; columnIndex < _rowSize; ++columnIndex)
				{
					std::copy(iter, (iter + rowIndex), newIter);
					iter += rowIndex;
					newIter += rowIndex;

					(*newIter) = (*rowVectorIter);
					++newIter;
					++rowVectorIter;

					std::copy(iter, (iter + remainedSize), newIter);
					iter += remainedSize;
					newIter += remainedSize;
				}
			}

			_valueArray.swap(newValueArray);
			_columnSize += 1;
		}
	}


	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "insertRow", "\"rowIndex\" is out of range." };
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::eraseColumn(const size_type columnIndex)
{
	if (isEmpty())
		throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "eraseColumn", "This is empty." };

	else
	{
		if (columnIndex < _rowSize)
		{
			std::valarray<value_type> newValueArray(_columnSize * (_rowSize - 1));
			{
				auto erasedIter = std::begin(_valueArray);
				auto newErasedIter = std::begin(newValueArray);
				{
					erasedIter += (_columnSize * columnIndex);
					newErasedIter += (_columnSize * columnIndex);
				}

				std::copy(std::begin(_valueArray), erasedIter, std::begin(newValueArray));
				std::copy((erasedIter + _columnSize), std::end(_valueArray), newErasedIter);
			}

			_valueArray.swap(newValueArray);
			--_rowSize;
		}


		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "eraseColumn", "\"columnIndex\" is out of range." };
	}
}

template <typename T>
inline void MathToolkit::LinearAlgebra::NumericalMatrix<T, 0, 0>::eraseRow(const size_type rowIndex)
{
	if (isEmpty())
		throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "eraseRow", "This is empty." };

	else
	{
		if (rowIndex < _columnSize)
		{
			std::valarray<value_type> newValueArray((_columnSize - 1) * _rowSize);
			{
				auto iter = std::begin(_valueArray);
				auto newIter = std::begin(newValueArray);
				size_type remainedSize = _columnSize - (rowIndex + 1);

				for (size_type columnIndex = 0; columnIndex < _rowSize; ++columnIndex)
				{
					std::copy(iter, (iter + rowIndex), newIter);
					iter += (rowIndex + 1);
					newIter += rowIndex;

					std::copy(iter, (iter + remainedSize), newIter);
					iter += remainedSize;
					newIter += remainedSize;
				}
			}

			_valueArray.swap(newValueArray);
			--_columnSize;
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "eraseRow", "\"rowIndex\" is out of range." };
	}
}

// Customize methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************

#endif // !MATHTOOLKIT_LINEARALGEBRA_NUMERICALMATRIX_H
