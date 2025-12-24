#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_COORDINATIONPOLYHEDRACONNECTIONPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_COORDINATIONPOLYHEDRACONNECTIONPARAMETERS_H

#include <cmath>

#include "ArgumentOutOfRangeException.h"

#include "StreamReader.h"

#include "GeometricalConstraintParameters.h"
#include "StructuralOptimizationParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			namespace Internal
			{
				class CoordinationPolyhedraConnectionParameters
				{
					using size_type = std::size_t;

					using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;
					using StructuralOptimizationParameters = MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters;


// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Constructors, destructor, and operators

				public:
					CoordinationPolyhedraConnectionParameters() noexcept;
					virtual ~CoordinationPolyhedraConnectionParameters() = default;

					CoordinationPolyhedraConnectionParameters(const CoordinationPolyhedraConnectionParameters&) = default;
					CoordinationPolyhedraConnectionParameters(CoordinationPolyhedraConnectionParameters&&) noexcept = default;
					CoordinationPolyhedraConnectionParameters& operator=(const CoordinationPolyhedraConnectionParameters&) = default;
					CoordinationPolyhedraConnectionParameters& operator=(CoordinationPolyhedraConnectionParameters&&) noexcept = default;

				// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Property

					const GeometricalConstraintParameters& geometricalConstraintParameters() const noexcept;
					const StructuralOptimizationParameters& molecularOptimizationParameters() const noexcept;

					void setGeometricalConstraintParameters(const GeometricalConstraintParameters&) noexcept;
					void setStructuralOptimizationParameters(const StructuralOptimizationParameters&) noexcept;

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
					void validateStructuralOptimizationParameters() const;

				// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

				private:
					GeometricalConstraintParameters _geometricalConstraintParameters;
					StructuralOptimizationParameters _molecularOptimizationParameters;


					static double s_defaultAttractiveForceConstant;
					static double s_defaultRepulsiveForceConstant;

					static size_type s_defaultMaxStructuralOptimizing;
					static double s_defaultInitialMaxAtomicDisplacement;
					static double s_defaultFinalMaxAtomicDisplacement;
					static double s_defaultFeasibleGeometricalConstraintErrorRate;
				};
			}
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

inline const MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::GeometricalConstraintParameters& MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::geometricalConstraintParameters() const noexcept
{
	return _geometricalConstraintParameters;
}

inline const MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::StructuralOptimizationParameters& MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::molecularOptimizationParameters() const noexcept
{
	return _molecularOptimizationParameters;
}

inline void MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::setGeometricalConstraintParameters(const GeometricalConstraintParameters& parameters) noexcept
{
	_geometricalConstraintParameters = parameters;
}

inline void MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::CoordinationPolyhedraConnectionParameters::setStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_molecularOptimizationParameters = parameters;
	validateStructuralOptimizationParameters();
}

// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_COORDINATIONPOLYHEDRACONNECTIONPARAMETERS_H
