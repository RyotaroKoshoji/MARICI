#ifndef MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYZER_H
#define MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYZER_H

#include "CrystalDesigner.h"

#include "CrystallographicStructure.h"

#include "OptimalCrystalStructure.h"
#include "ConstrainingCrystalStructure.h"

#include "CrystalOptimalityAnalysisParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace Analysis
	{
		class CrystalOptimalityAnalyzer
		{
			using size_type = std::size_t;
			using ConstrainingCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure;

			using CrystalDesigner = MathematicalCrystalChemistry::Design::CrystalDesigner;

// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalOptimalityAnalyzer() noexcept;
			explicit CrystalOptimalityAnalyzer(const CrystalOptimalityAnalysisParameters&);

			virtual ~CrystalOptimalityAnalyzer() = default;

			CrystalOptimalityAnalyzer(const CrystalOptimalityAnalyzer&) = default;
			CrystalOptimalityAnalyzer(CrystalOptimalityAnalyzer&&) noexcept = default;
			CrystalOptimalityAnalyzer& operator=(const CrystalOptimalityAnalyzer&) = default;
			CrystalOptimalityAnalyzer& operator=(CrystalOptimalityAnalyzer&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			bool isActive() const noexcept;

			void setOptimalityAnalysisParameters(const CrystalOptimalityAnalysisParameters&) noexcept;

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			bool isFeasible(ConstrainingCrystalStructure&) const;
			void optimizeStructure(ConstrainingCrystalStructure&) const;

			template <typename A>
			bool isFeasible(const ChemToolkit::Crystallography::CrystallographicStructure<A>&) const;

			template <typename A>
			ConstrainingCrystalStructure optimizeStructure(const ChemToolkit::Crystallography::CrystallographicStructure<A>&) const;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			bool _isActive;
			CrystalDesigner _crystalDesigner;
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

inline bool MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::isActive() const noexcept
{
	return _isActive;
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::setOptimalityAnalysisParameters(const CrystalOptimalityAnalysisParameters& parameters) noexcept
{
	_isActive = parameters.needAnalysis();

	MathematicalCrystalChemistry::Design::CrystalDesignParameters crystalDesignParameters;
	{
		crystalDesignParameters.setGeometricalConstraintParameters(parameters.geometricalConstraintParameters());
		crystalDesignParameters.setGlobalStructuralOptimizationParameters(parameters.globalStructuralOptimizationParameters());
		crystalDesignParameters.setLocalStructuralOptimizationParameters(parameters.localStructuralOptimizationParameters());
		crystalDesignParameters.setPreciseStructuralOptimizationParameters(parameters.preciseStructuralOptimizationParameters());


		size_type maxTotalStructuralOptimizing = 0;
		{
			maxTotalStructuralOptimizing += parameters.globalStructuralOptimizationParameters().maxStructuralOptimizing();
			maxTotalStructuralOptimizing += parameters.localStructuralOptimizationParameters().maxStructuralOptimizing();
			maxTotalStructuralOptimizing += parameters.preciseStructuralOptimizationParameters().maxStructuralOptimizing();
		}
		crystalDesignParameters.setMaxTotalStructuralOptimizing(maxTotalStructuralOptimizing);
	}

	_crystalDesigner.setCrystalDesignParameters(crystalDesignParameters);
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

inline bool MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::isFeasible(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	optimizeStructure(constrainingCrystalStructure);
	constrainingCrystalStructure.setFeasibleErrorRate(_crystalDesigner.preciseStructuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());

	return constrainingCrystalStructure.isFeasible();
}

inline void MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::optimizeStructure(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	_crystalDesigner.execute(constrainingCrystalStructure);
}

template <typename A>
inline bool MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::isFeasible(const ChemToolkit::Crystallography::CrystallographicStructure<A>& crystalStructure) const
{
	ConstrainingCrystalStructure constrainingCrystalStructure;
	{
		using MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom;
		using MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure;

		std::vector<OptimalAtom> optimalAtoms;
		{
			for (const auto& atom : crystalStructure.atoms())
				optimalAtoms.push_back(OptimalAtom{ atom.ionicAtomicNumber(), atom.cartesianCoordinate() });
		}


		OptimalCrystalStructure optimalCrystalStructure{ crystalStructure.unitCell(), std::move(optimalAtoms) };
		constrainingCrystalStructure = ConstrainingCrystalStructure{ optimalCrystalStructure };
	}

	return isFeasible(constrainingCrystalStructure);
}

template <typename A>
inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure MathematicalCrystalChemistry::Analysis::CrystalOptimalityAnalyzer::optimizeStructure(const ChemToolkit::Crystallography::CrystallographicStructure<A>& structure) const
{
	using MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure;
	OptimalCrystalStructure optimalCrystalStructure;
	{
		using MathematicalCrystalChemistry::CrystalModel::Components::OptimalAtom;
		std::vector<OptimalAtom> optimalAtoms;
		{
			for (const auto& atom : structure.atoms())
				optimalAtoms.push_back(OptimalAtom{ atom.ionicAtomicNumber(), atom.cartesianCoordinate() });
		}

		optimalCrystalStructure = OptimalCrystalStructure{ structure.unitCell(), std::move(optimalAtoms) };
	}

	ConstrainingCrystalStructure constrainingCrystalStructure{ optimalCrystalStructure };
	optimizeStructure(constrainingCrystalStructure);


	return constrainingCrystalStructure;
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_ANALYSIS_CRYSTALOPTIMALITYANALYZER_H
