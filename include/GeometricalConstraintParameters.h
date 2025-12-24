#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_GEOMETRICALCONSTRAINTPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_GEOMETRICALCONSTRAINTPARAMETERS_H

#include "ArgumentOutOfRangeException.h"

#include "StreamReader.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Optimization
		{
			class GeometricalConstraintParameters
			{
				using size_type = std::size_t;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				GeometricalConstraintParameters() noexcept;
				virtual ~GeometricalConstraintParameters() = default;

				GeometricalConstraintParameters(const GeometricalConstraintParameters&) = default;
				GeometricalConstraintParameters(GeometricalConstraintParameters&&) noexcept = default;
				GeometricalConstraintParameters& operator=(const GeometricalConstraintParameters&) = default;
				GeometricalConstraintParameters& operator=(GeometricalConstraintParameters&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				size_type interatomicDistanceTracerTimeout() const noexcept;
				size_type unitCellReductionTimeout() const noexcept;
				double minimumExclusionDistanceRatio() const noexcept;
				double interatomicDistanceTracerCutoffRatio() const noexcept;
				double interatomicDistanceConstrainerCutoffRatio() const noexcept;

				static size_type defaultInteratomicDistanceTracerTimeout() noexcept;
				static size_type defaultUnitCellReductionTimeout() noexcept;
				static double defaultMinimumExclusionDistanceRatio() noexcept;
				static double defaultInteratomicDistanceTracerCutoffRatio() noexcept;
				static double defaultInteratomicDistanceConstrainerCutoffRatio() noexcept;


				void setInteratomicDistanceTracerTimeout() noexcept;
				void setInteratomicDistanceTracerTimeout(const size_type);
				void setUnitCellReductionTimeout() noexcept;
				void setUnitCellReductionTimeout(const size_type);
				void setMinimumExclusionDistanceRatio() noexcept;
				void setMinimumExclusionDistanceRatio(const double);
				void setInteratomicDistanceTracerCutoffRatio() noexcept;
				void setInteratomicDistanceTracerCutoffRatio(const double);
				void setInteratomicDistanceConstrainerCutoffRatio() noexcept;
				void setInteratomicDistanceConstrainerCutoffRatio(const double);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void initialize() noexcept;
				void initialize(const System::IO::StreamReader& streamReader);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void validateInitializedValues();

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				size_type _interatomicDistanceTracerTimeout;
				size_type _unitCellReductionTimeout;
				double _minimumExclusionDistanceRatio;
				double _interatomicDistanceTracerCutoffRatio;
				double _interatomicDistanceConstrainerCutoffRatio;

				static size_type s_defaultInteratomicDistanceTracerTimeout;
				static size_type s_defaultUnitCellReductionTimeout;
				static double s_defaultMinimumExclusionDistanceRatio;
				static double s_defaultInteratomicDistanceTracerCutoffRatio;
				static double s_defaultInteratomicDistanceConstrainerCutoffRatio;
			};
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
// Property

inline MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::size_type MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::interatomicDistanceTracerTimeout() const noexcept
{
	return _interatomicDistanceTracerTimeout;
}

inline MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::size_type MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::unitCellReductionTimeout() const noexcept
{
	return _unitCellReductionTimeout;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::minimumExclusionDistanceRatio() const noexcept
{
	return _minimumExclusionDistanceRatio;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::interatomicDistanceTracerCutoffRatio() const noexcept
{
	return _interatomicDistanceTracerCutoffRatio;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::interatomicDistanceConstrainerCutoffRatio() const noexcept
{
	return _interatomicDistanceConstrainerCutoffRatio;
}

inline MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::size_type MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerTimeout() noexcept
{
	return s_defaultInteratomicDistanceTracerTimeout;
}

inline MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::size_type MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultUnitCellReductionTimeout() noexcept
{
	return s_defaultUnitCellReductionTimeout;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultMinimumExclusionDistanceRatio() noexcept
{
	return s_defaultMinimumExclusionDistanceRatio;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerCutoffRatio() noexcept
{
	return s_defaultInteratomicDistanceTracerCutoffRatio;
}

inline double MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceConstrainerCutoffRatio() noexcept
{
	return s_defaultInteratomicDistanceConstrainerCutoffRatio;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceTracerTimeout() noexcept
{
	_interatomicDistanceTracerTimeout = s_defaultInteratomicDistanceTracerTimeout;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceTracerTimeout(const size_type val)
{
	if (0 < val)
		_interatomicDistanceTracerTimeout = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInteratomicDistanceTracerTimeout", "Input value is zero." };
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setUnitCellReductionTimeout() noexcept
{
	_unitCellReductionTimeout = s_defaultUnitCellReductionTimeout;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setUnitCellReductionTimeout(const size_type val)
{
	if (0 < val)
		_unitCellReductionTimeout = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setUnitCellReductionTimeout", "Input value is zero." };
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setMinimumExclusionDistanceRatio() noexcept
{
	_minimumExclusionDistanceRatio = s_defaultMinimumExclusionDistanceRatio;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setMinimumExclusionDistanceRatio(const double val)
{
	if (1.0 < val)
		_minimumExclusionDistanceRatio = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMinimumExclusionDistanceRatio", "Input value is not more than zero." };
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceTracerCutoffRatio() noexcept
{
	_interatomicDistanceTracerCutoffRatio = s_defaultInteratomicDistanceTracerCutoffRatio;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceTracerCutoffRatio(const double val)
{
	if (0.0 < val)
		_interatomicDistanceTracerCutoffRatio = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInteratomicDistanceTrackerCutoffRatio", "Input value is not more than zero." };
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceConstrainerCutoffRatio() noexcept
{
	_interatomicDistanceConstrainerCutoffRatio = s_defaultInteratomicDistanceConstrainerCutoffRatio;
}

inline void MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::setInteratomicDistanceConstrainerCutoffRatio(const double val)
{
	if (0.0 < val)
		_interatomicDistanceConstrainerCutoffRatio = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInteratomicDistanceConstrainerCutoffRatio", "Input value is not more than zero." };
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_GEOMETRICALCONSTRAINTPARAMETERS_H
