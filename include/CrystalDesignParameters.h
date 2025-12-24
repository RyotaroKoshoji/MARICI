#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNPARAMETERS_H

#include <vector>
#include <filesystem>

#include "ArgumentOutOfRangeException.h"

#include "StreamReader.h"

#include "InitialStructureGenerationParameters.h"
#include "GeometricalConstraintParameters.h"
#include "StructuralOptimizationParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		class CrystalDesignParameters
		{
			using size_type = std::size_t;

			using InitialStructureGenerationParameters = MathematicalCrystalChemistry::Design::Generation::InitialStructureGenerationParameters;
			using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;
			using StructuralOptimizationParameters = MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalDesignParameters() noexcept;
			virtual ~CrystalDesignParameters() = default;

			CrystalDesignParameters(const CrystalDesignParameters&) = default;
			CrystalDesignParameters(CrystalDesignParameters&&) noexcept = default;
			CrystalDesignParameters& operator=(const CrystalDesignParameters&) = default;
			CrystalDesignParameters& operator=(CrystalDesignParameters&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			size_type maxTotalStructuralOptimizing() const noexcept;
			size_type maxCeaselessGlobalStructuralOptimizing() const noexcept;
			const InitialStructureGenerationParameters& initialStructureGenerationParameters() const noexcept;
			const GeometricalConstraintParameters& geometricalConstraintParameters() const noexcept;

			const StructuralOptimizationParameters& globalStructuralOptimizationParameters() const noexcept;
			const StructuralOptimizationParameters& localStructuralOptimizationParameters() const noexcept;
			const StructuralOptimizationParameters& preciseStructuralOptimizationParameters() const noexcept;


			void setMaxTotalStructuralOptimizing(const size_type) noexcept;
			void setMaxCeaselessGlobalStructuralOptimizing() noexcept;
			void setMaxCeaselessGlobalStructuralOptimizing(const size_type);
			void setInitialStructureGenerationParameters(const InitialStructureGenerationParameters&) noexcept;
			void setGeometricalConstraintParameters(const GeometricalConstraintParameters&) noexcept;

			void setGlobalStructuralOptimizationParameters(const StructuralOptimizationParameters&) noexcept;
			void setLocalStructuralOptimizationParameters(const StructuralOptimizationParameters&) noexcept;
			void setPreciseStructuralOptimizationParameters(const StructuralOptimizationParameters&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			void initialize() noexcept;
			void initialize(const System::IO::StreamReader&);

			bool isValid() const noexcept;

			// Methods
	// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			size_type _maxTotalStructuralOptimizing;
			size_type _maxCeaselessGlobalStructuralOptimizing;

			InitialStructureGenerationParameters _initialStructureGenerationParameters;
			GeometricalConstraintParameters _geometricalConstraintParameters;

			StructuralOptimizationParameters _globalStructuralOptimizationParameters;
			StructuralOptimizationParameters _localStructuralOptimizationParameters;
			StructuralOptimizationParameters _preciseStructuralOptimizationParameters;


			static size_type s_defaultMaxCeaselessGlobalStructuralOptimizing;
		};
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

inline MathematicalCrystalChemistry::Design::CrystalDesignParameters::size_type MathematicalCrystalChemistry::Design::CrystalDesignParameters::maxTotalStructuralOptimizing() const noexcept
{
	return _maxTotalStructuralOptimizing;
}

inline MathematicalCrystalChemistry::Design::CrystalDesignParameters::size_type MathematicalCrystalChemistry::Design::CrystalDesignParameters::maxCeaselessGlobalStructuralOptimizing() const noexcept
{
	return _maxCeaselessGlobalStructuralOptimizing;
}

inline const MathematicalCrystalChemistry::Design::Generation::InitialStructureGenerationParameters& MathematicalCrystalChemistry::Design::CrystalDesignParameters::initialStructureGenerationParameters() const noexcept
{
	return _initialStructureGenerationParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters& MathematicalCrystalChemistry::Design::CrystalDesignParameters::geometricalConstraintParameters() const noexcept
{
	return _geometricalConstraintParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::CrystalDesignParameters::globalStructuralOptimizationParameters() const noexcept
{
	return _globalStructuralOptimizationParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::CrystalDesignParameters::localStructuralOptimizationParameters() const noexcept
{
	return _localStructuralOptimizationParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::CrystalDesignParameters::preciseStructuralOptimizationParameters() const noexcept
{
	return _preciseStructuralOptimizationParameters;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setMaxTotalStructuralOptimizing(const size_type val) noexcept
{
	_maxTotalStructuralOptimizing = val;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setMaxCeaselessGlobalStructuralOptimizing() noexcept
{
	_maxCeaselessGlobalStructuralOptimizing = s_defaultMaxCeaselessGlobalStructuralOptimizing;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setMaxCeaselessGlobalStructuralOptimizing(const size_type val)
{
	if (0 < val)
		_maxCeaselessGlobalStructuralOptimizing = val;
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "setMaxCeaselessGlobalStructuralOptimizing", "Argument value is zero." };
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setInitialStructureGenerationParameters(const InitialStructureGenerationParameters& parameters) noexcept
{
	_initialStructureGenerationParameters = parameters;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setGeometricalConstraintParameters(const GeometricalConstraintParameters& parameters) noexcept
{
	_geometricalConstraintParameters = parameters;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setGlobalStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_globalStructuralOptimizationParameters = parameters;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setLocalStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_localStructuralOptimizationParameters = parameters;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesignParameters::setPreciseStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_preciseStructuralOptimizationParameters = parameters;
}

// Property
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
// Methods

inline bool MathematicalCrystalChemistry::Design::CrystalDesignParameters::isValid() const noexcept
{
	if (!(_initialStructureGenerationParameters.isValid()))
		return false;

	return true;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNPARAMETERS_H
