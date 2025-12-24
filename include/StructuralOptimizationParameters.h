#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_STRUCTURALOPTIMIZATIONPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_STRUCTURALOPTIMIZATIONPARAMETERS_H

#include <cmath>

#include "GeometricalConstraintParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		namespace Optimization
		{
			class StructuralOptimizationParameters
			{
				using size_type = std::size_t;

			public:
				enum class OptimizationType
				{
					global,
					local,
					precise
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				StructuralOptimizationParameters() noexcept;
				virtual ~StructuralOptimizationParameters() = default;

				StructuralOptimizationParameters(const StructuralOptimizationParameters&) = default;
				StructuralOptimizationParameters(StructuralOptimizationParameters&&) noexcept = default;
				StructuralOptimizationParameters& operator=(const StructuralOptimizationParameters&) = default;
				StructuralOptimizationParameters& operator=(StructuralOptimizationParameters&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				double pressure() const noexcept;
				double attractiveForceConstant() const noexcept;
				double repulsiveForceConstant() const noexcept;

				size_type maxStructuralOptimizing() const noexcept;
				double initialMaxAtomicDisplacement() const noexcept;
				double initialMaxUnitCellDisplacement() const noexcept;
				double displacementDecreasingFactor() const noexcept;
				double feasibleGeometricalConstraintErrorRate() const noexcept;


				void setPressure(const double);
				void setAttractiveForceConstant(const double);
				void setRepulsiveForceConstant(const double);

				void setMaxStructuralOptimizing(const size_type);
				void setMaxAtomicDisplacements(const double initialMaxAtomicDisplacement, const double finalMaxAtomicDisplacement);
				void setInitialMaxUnitCellDisplacement(const double maxUnitCellDisplacementFactor);
				void setFeasibleGeometricalConstraintErrorRate(const double);

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void initialize() noexcept;
				void initialize(const OptimizationType);
				void initialize(const System::IO::StreamReader& streamReader, const OptimizationType);

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				void validateInitializedValues() const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				double _pressure;
				double _attractiveForceConstant;
				double _repulsiveForceConstant;

				size_type _maxStructuralOptimizing;
				double _initialMaxAtomicDisplacement;
				double _initialMaxUnitCellDisplacement;
				double _displacementDecreasingFactor;
				double _feasibleGeometricalConstraintErrorRate;


				static double s_defaultPressure;
				static double s_defaultAttractiveForceConstant;
				static double s_defaultRepulsiveForceConstant;
				static double s_defaultMaxUnitCellDisplacementFactor;

				static size_type s_defaultGlobalMaxStructuralOptimizing;
				static double s_defaultGlobalInitialMaxAtomicDisplacement;
				static double s_defaultGlobalFinalMaxAtomicDisplacement;
				static double s_defaultGlobalFeasibleGeometricalConstraintErrorRate;

				static size_type s_defaultLocalMaxStructuralOptimizing;
				static double s_defaultLocalInitialMaxAtomicDisplacement;
				static double s_defaultLocalFinalMaxAtomicDisplacement;
				static double s_defaultLocalFeasibleGeometricalConstraintErrorRate;

				static size_type s_defaultPreciseMaxStructuralOptimizing;
				static double s_defaultPreciseInitialMaxAtomicDisplacement;
				static double s_defaultPreciseFinalMaxAtomicDisplacement;
				static double s_defaultPreciseFeasibleGeometricalConstraintErrorRate;
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

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::pressure() const noexcept
{
	return _pressure;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::attractiveForceConstant() const noexcept
{
	return _attractiveForceConstant;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::repulsiveForceConstant() const noexcept
{
	return _repulsiveForceConstant;
}

inline MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::size_type MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::maxStructuralOptimizing() const noexcept
{
	return _maxStructuralOptimizing;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::initialMaxAtomicDisplacement() const noexcept
{
	return _initialMaxAtomicDisplacement;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::initialMaxUnitCellDisplacement() const noexcept
{
	return _initialMaxUnitCellDisplacement;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::displacementDecreasingFactor() const noexcept
{
	return _displacementDecreasingFactor;
}

inline double MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::feasibleGeometricalConstraintErrorRate() const noexcept
{
	return _feasibleGeometricalConstraintErrorRate;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setPressure(const double val)
{
	if (val < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setPressure", "Argument value is less than zero." };
	else
		_pressure = val;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setAttractiveForceConstant(const double val)
{
	if (val < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setAttractiveForceConstant", "Argument value is less than zero." };
	else
		_attractiveForceConstant = val;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setRepulsiveForceConstant(const double val)
{
	if (0.0 < val)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setRepulsiveForceConstant", "Argument value is more than zero." };
	else
		_repulsiveForceConstant = val;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setMaxStructuralOptimizing(const size_type value)
{
	double finalDisplacement = _initialMaxAtomicDisplacement * std::pow(_displacementDecreasingFactor, static_cast<double>(_maxStructuralOptimizing));
	_maxStructuralOptimizing = value;


	if ((0 < _maxStructuralOptimizing) && (0.0 < _initialMaxAtomicDisplacement))
		_displacementDecreasingFactor = std::pow((finalDisplacement / _initialMaxAtomicDisplacement), (1.0 / static_cast<double>(_maxStructuralOptimizing)));
	else
		_displacementDecreasingFactor = 0.0;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setMaxAtomicDisplacements(const double initialValue, const double finalValue)
{
	if (initialValue < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMaxAtomicDisplacements", "Initial maximum displacement is less than zero." };

	else
	{
		if (finalValue < 0.0)
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMaxAtomicDisplacements", "Final maximum displacement is less than zero." };

		else
		{
			if (initialValue < finalValue)
				throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMaxAtomicDisplacements", "Initial maximum displacement is less than the final maximum displacement." };

			else
			{
				double maxUnitCellDisplacementFactor = 0.0;
				{
					if (0.0 < _initialMaxAtomicDisplacement)
						maxUnitCellDisplacementFactor = (_initialMaxUnitCellDisplacement / _initialMaxAtomicDisplacement);
				}

				_initialMaxAtomicDisplacement = initialValue;
				_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * maxUnitCellDisplacementFactor;


				if ((0 < _maxStructuralOptimizing) && (0.0 < _initialMaxAtomicDisplacement))
					_displacementDecreasingFactor = std::pow((finalValue / _initialMaxAtomicDisplacement), (1.0 / static_cast<double>(_maxStructuralOptimizing)));
				else
					_displacementDecreasingFactor = 0.0;
			}
		}
	}
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setInitialMaxUnitCellDisplacement(const double maxUnitCellDisplacementFactor)
{
	if (maxUnitCellDisplacementFactor < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setInitialMaxUnitCellDisplacement", "Maximum unit cell displacement factor is less than zero." };
	else
		_initialMaxUnitCellDisplacement = _initialMaxAtomicDisplacement * maxUnitCellDisplacementFactor;
}

inline void MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters::setFeasibleGeometricalConstraintErrorRate(const double val)
{
	if (val < 0.0)
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setFeasibleGeometricalConstraintErrorRate", "Geometrical constraint error rate is less than zero." };
	else
		_feasibleGeometricalConstraintErrorRate = val;
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_OPTIMIZATION_STRUCTURALOPTIMIZATIONPARAMETERS_H
