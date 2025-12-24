#ifndef MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYSISPARAMETERS_H
#define MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYSISPARAMETERS_H

#include "StreamReader.h"

#include "GeometricalConstraintParameters.h"
#include "StructuralOptimizationParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Analysis
	{
		class CrystalOptimalityAnalysisParameters
		{
			using size_type = std::size_t;
			using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;
			using StructuralOptimizationParameters = MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalOptimalityAnalysisParameters() noexcept;
			virtual ~CrystalOptimalityAnalysisParameters() = default;

			CrystalOptimalityAnalysisParameters(const CrystalOptimalityAnalysisParameters&) = default;
			CrystalOptimalityAnalysisParameters(CrystalOptimalityAnalysisParameters&&) noexcept = default;
			CrystalOptimalityAnalysisParameters& operator=(const CrystalOptimalityAnalysisParameters&) = default;
			CrystalOptimalityAnalysisParameters& operator=(CrystalOptimalityAnalysisParameters&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			bool needAnalysis() const noexcept;

			const GeometricalConstraintParameters& geometricalConstraintParameters() const noexcept;
			const StructuralOptimizationParameters& globalStructuralOptimizationParameters() const noexcept;
			const StructuralOptimizationParameters& localStructuralOptimizationParameters() const noexcept;
			const StructuralOptimizationParameters& preciseStructuralOptimizationParameters() const noexcept;


			void setAnalysisNecessity(const bool) noexcept;

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

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Private methods

		private:
			bool toNecessity(const std::string&) const;

		// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			bool _needAnalysis;

			GeometricalConstraintParameters _geometricalConstraintParameters;
			StructuralOptimizationParameters _globalStructuralOptimizationParameters;
			StructuralOptimizationParameters _localStructuralOptimizationParameters;
			StructuralOptimizationParameters _preciseStructuralOptimizationParameters;
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

inline bool MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::needAnalysis() const noexcept
{
	return _needAnalysis;
}

inline const MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters& MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::geometricalConstraintParameters() const noexcept
{
	return _geometricalConstraintParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::globalStructuralOptimizationParameters() const noexcept
{
	return _globalStructuralOptimizationParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::localStructuralOptimizationParameters() const noexcept
{
	return _localStructuralOptimizationParameters;
}

inline const MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters& MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::preciseStructuralOptimizationParameters() const noexcept
{
	return _preciseStructuralOptimizationParameters;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::setAnalysisNecessity(const bool necessity) noexcept
{
	_needAnalysis = necessity;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::setGeometricalConstraintParameters(const GeometricalConstraintParameters& parameters) noexcept
{
	_geometricalConstraintParameters = parameters;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::setGlobalStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_globalStructuralOptimizationParameters = parameters;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::setLocalStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
{
	_localStructuralOptimizationParameters = parameters;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalysisParameters::setPreciseStructuralOptimizationParameters(const StructuralOptimizationParameters& parameters) noexcept
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


#endif // !MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYSISPARAMETERS_H
