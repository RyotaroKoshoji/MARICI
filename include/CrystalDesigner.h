#ifndef MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNER_H
#define MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNER_H

#include <utility>

#include "ChemicalComposition.h"

#include "CrystalDesignParameters.h"

#include "RandomStructureGenerator.h"
#include "CrystalOptimizer.h"
#include "CrystalDesignRecorder.h"

#include "ConstrainingAtomicSpecies.h"
#include "ConstrainingCrystalStructure.h"
#include "ObjectiveCrystalStructure.h"


namespace MathematicalCrystalChemistry
{
	namespace Design
	{
		class CrystalDesigner
		{
			using size_type = std::size_t;

			using GeometricalConstraintParameters = MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters;
			using StructuralOptimizationParameters = MathematicalCrystalChemistry::Design::Optimization::StructuralOptimizationParameters;
			using RandomStructureGenerator = MathematicalCrystalChemistry::Design::Generation::RandomStructureGenerator;
			using CrystalOptimizer = MathematicalCrystalChemistry::Design::Optimization::CrystalOptimizer;
			using CrystalDesignRecorder = MathematicalCrystalChemistry::Design::Diagnostics::CrystalDesignRecorder;

			using ConstrainingAtomicSpecies = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies;
			using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<ConstrainingAtomicSpecies>;			
			using ConstrainingCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure;
			using ObjectiveCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveCrystalStructure;


// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Constructors, destructor, and operators

		public:
			CrystalDesigner() noexcept;
			explicit CrystalDesigner(const CrystalDesignParameters&);

			virtual ~CrystalDesigner() = default;

			CrystalDesigner(const CrystalDesigner&) = default;
			CrystalDesigner(CrystalDesigner&&) noexcept = default;
			CrystalDesigner& operator=(const CrystalDesigner&) = default;
			CrystalDesigner& operator=(CrystalDesigner&&) noexcept = default;

		// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Property

			const StructuralOptimizationParameters& preciseStructuralOptimizationParameters() const noexcept;

			void setCrystalDesignParameters(const CrystalDesignParameters&);
			void setProductChemicalComposition();
			void setProductChemicalComposition(const ChemicalComposition&);

		// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Methods

			ConstrainingCrystalStructure execute() const;
			ConstrainingCrystalStructure execute(CrystalDesignRecorder&) const;

			void execute(ConstrainingCrystalStructure& initialStructure) const;
			void execute(ConstrainingCrystalStructure& initialStructure, CrystalDesignRecorder&) const;

		// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
		// Private methods

		private:
			void initializeTimers() const noexcept;

			void applyGlobalStructuralOptimization(ObjectiveCrystalStructure&) const;
			void applyGlobalStructuralOptimization(ObjectiveCrystalStructure&, CrystalDesignRecorder&) const;

			bool applyLocalStructuralOptimization(ObjectiveCrystalStructure&) const;
			bool applyLocalStructuralOptimization(ObjectiveCrystalStructure&, CrystalDesignRecorder&) const;

			bool applyPreciseStructuralOptimization(ObjectiveCrystalStructure&) const;
			bool applyPreciseStructuralOptimization(ObjectiveCrystalStructure&, CrystalDesignRecorder&) const;

			bool isFeasible(ConstrainingCrystalStructure&) const;

			void reduceStructure(ConstrainingCrystalStructure&) const;
			void updateConstraints(ConstrainingCrystalStructure&) const;
			void forceUpdateConstraints(ConstrainingCrystalStructure&) const;

		// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

		private:
			size_type _maxTotalStructuralOptimizing;
			size_type _maxCeaselessGlobalStructuralOptimizing;
			GeometricalConstraintParameters _geometricalConstraintParameters;

			RandomStructureGenerator _randomStructureGenerator;
			CrystalOptimizer _globalStructuralOptimizer;
			CrystalOptimizer _localStructuralOptimizer;
			CrystalOptimizer _preciseStructuralOptimizer;


			mutable size_type m_totalStructuralOptimizing;
			mutable size_type m_ceaselessGlobalStructuralOptimizing;
			mutable size_type m_interatomicDistanceTrackerUsing;
			mutable size_type m_unitCellUsing;
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

inline const MathematicalCrystalChemistry::Design::CrystalDesigner::StructuralOptimizationParameters& MathematicalCrystalChemistry::Design::CrystalDesigner::preciseStructuralOptimizationParameters() const noexcept
{
	return _preciseStructuralOptimizer.structuralOptimizationParameters();
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::setCrystalDesignParameters(const CrystalDesignParameters& parameters)
{
	_maxTotalStructuralOptimizing = parameters.maxTotalStructuralOptimizing();
	_maxCeaselessGlobalStructuralOptimizing = parameters.maxCeaselessGlobalStructuralOptimizing();
	_geometricalConstraintParameters = parameters.geometricalConstraintParameters();

	_randomStructureGenerator.setRandomStructureGenerationParameters(parameters.initialStructureGenerationParameters().randomStructureGenerationParameters());
	_globalStructuralOptimizer.setParameters(parameters.globalStructuralOptimizationParameters(), parameters.geometricalConstraintParameters());
	_localStructuralOptimizer.setParameters(parameters.localStructuralOptimizationParameters(), parameters.geometricalConstraintParameters());
	_preciseStructuralOptimizer.setParameters(parameters.preciseStructuralOptimizationParameters(), parameters.geometricalConstraintParameters());
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::setProductChemicalComposition()
{
	_randomStructureGenerator.setGeneratingChemicalComposition();
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::setProductChemicalComposition(const ChemicalComposition& composition)
{
	_randomStructureGenerator.setGeneratingChemicalComposition(composition);
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

inline MathematicalCrystalChemistry::Design::CrystalDesigner::ConstrainingCrystalStructure MathematicalCrystalChemistry::Design::CrystalDesigner::execute() const
{
	ConstrainingCrystalStructure constrainingCrystalStructure = _randomStructureGenerator.next();
	execute(constrainingCrystalStructure);

	return constrainingCrystalStructure;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure MathematicalCrystalChemistry::Design::CrystalDesigner::execute(CrystalDesignRecorder& crystalDesignRecorder) const
{
	ConstrainingCrystalStructure constrainingCrystalStructure = _randomStructureGenerator.next();
	execute(constrainingCrystalStructure, crystalDesignRecorder);

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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private methods

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::initializeTimers() const noexcept
{
	m_totalStructuralOptimizing = 0;
	m_ceaselessGlobalStructuralOptimizing = 0;
	m_interatomicDistanceTrackerUsing = 0;
	m_unitCellUsing = 0;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::applyGlobalStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure) const
{
	_globalStructuralOptimizer.execute(objectiveCrystalStructure);

	size_type structuralOptimizing = _globalStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;
	m_ceaselessGlobalStructuralOptimizing += structuralOptimizing;
	m_interatomicDistanceTrackerUsing += structuralOptimizing;
	m_unitCellUsing += structuralOptimizing;
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::applyGlobalStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure, CrystalDesignRecorder& recorder) const
{
	_globalStructuralOptimizer.execute(objectiveCrystalStructure, recorder);

	size_type structuralOptimizing = _globalStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;
	m_ceaselessGlobalStructuralOptimizing += structuralOptimizing;
	m_interatomicDistanceTrackerUsing += structuralOptimizing;
	m_unitCellUsing += structuralOptimizing;
}

inline bool MathematicalCrystalChemistry::Design::CrystalDesigner::applyLocalStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure) const
{
	const double localErrorRate = _localStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate();
	_localStructuralOptimizer.execute(objectiveCrystalStructure);

	size_type structuralOptimizing = _localStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;


	if (objectiveCrystalStructure.isFeasible(localErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio()))
		return true;

	else
	{
		_localStructuralOptimizer.execute(objectiveCrystalStructure);
		m_totalStructuralOptimizing += structuralOptimizing;

		return objectiveCrystalStructure.isFeasible(localErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio());
	}
}

inline bool MathematicalCrystalChemistry::Design::CrystalDesigner::applyLocalStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure, CrystalDesignRecorder& recorder) const
{
	const double localErrorRate = _localStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate();
	_localStructuralOptimizer.execute(objectiveCrystalStructure, recorder);

	size_type structuralOptimizing = _localStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;


	if (objectiveCrystalStructure.isFeasible(localErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio()))
		return true;

	else
	{
		_localStructuralOptimizer.execute(objectiveCrystalStructure, recorder);
		m_totalStructuralOptimizing += structuralOptimizing;

		return objectiveCrystalStructure.isFeasible(localErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio());
	}
}

inline bool MathematicalCrystalChemistry::Design::CrystalDesigner::applyPreciseStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure) const
{
	const double preciseErrorRate = _preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate();
	_preciseStructuralOptimizer.execute(objectiveCrystalStructure);

	size_type structuralOptimizing = _preciseStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;


	return objectiveCrystalStructure.isFeasible(preciseErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio());
}

inline bool MathematicalCrystalChemistry::Design::CrystalDesigner::applyPreciseStructuralOptimization(ObjectiveCrystalStructure& objectiveCrystalStructure, CrystalDesignRecorder& recorder) const
{
	const double preciseErrorRate = _preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate();
	_preciseStructuralOptimizer.execute(objectiveCrystalStructure, recorder);

	size_type structuralOptimizing = _preciseStructuralOptimizer.structuralOptimizationParameters().maxStructuralOptimizing();
	m_totalStructuralOptimizing += structuralOptimizing;


	return objectiveCrystalStructure.isFeasible(preciseErrorRate, _geometricalConstraintParameters.minimumExclusionDistanceRatio());
}

inline bool MathematicalCrystalChemistry::Design::CrystalDesigner::isFeasible(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	constrainingCrystalStructure.reduceStructure();
	{
		if (constrainingCrystalStructure.hasFeasibleUnitCell())
		{
			constrainingCrystalStructure.updateTracingIndexPairs();
			constrainingCrystalStructure.createInteratomicDistanceConstraints();
			constrainingCrystalStructure.eraseInfeasibleChemicalBonds();

			m_interatomicDistanceTrackerUsing = 0;
			m_unitCellUsing = 0;

			return constrainingCrystalStructure.isFeasible();
		}

		else
			throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "forceUpdateConstraints", "Unit cell is infeasible despite the reduction." };
	}
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::reduceStructure(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	constrainingCrystalStructure.reduceStructure();
	{
		if (constrainingCrystalStructure.hasFeasibleUnitCell())
		{
			constrainingCrystalStructure.updateTracingIndexPairs();
			constrainingCrystalStructure.createInteratomicDistanceConstraints();
			constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();

			m_interatomicDistanceTrackerUsing = 0;
			m_unitCellUsing = 0;
		}

		else
			throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "reduceStructure", "Unit cell is infeasible despite the reduction." };
	}
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::updateConstraints(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	if (_geometricalConstraintParameters.unitCellReductionTimeout() < m_unitCellUsing)
	{
		constrainingCrystalStructure.reduceStructure();
		{
			if (constrainingCrystalStructure.hasFeasibleUnitCell())
			{
				constrainingCrystalStructure.updateTracingIndexPairs();
				constrainingCrystalStructure.createInteratomicDistanceConstraints();
				constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();

				m_interatomicDistanceTrackerUsing = 0;
				m_unitCellUsing = 0;
			}

			else
				throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "reduceStructure", "Unit cell is infeasible despite the reduction." };
		}
	}

	else
	{
		if (constrainingCrystalStructure.hasFeasibleUnitCell())
		{
			if (_geometricalConstraintParameters.interatomicDistanceTracerTimeout() < m_interatomicDistanceTrackerUsing)
			{
				constrainingCrystalStructure.normalizeFractionalCoordinates();
				constrainingCrystalStructure.updateTracingIndexPairs();
				constrainingCrystalStructure.createInteratomicDistanceConstraints();
				constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();

				m_interatomicDistanceTrackerUsing = 0;
			}

			else
			{
				constrainingCrystalStructure.createInteratomicDistanceConstraints();
				constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();
			}
		}

		else
			throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "updateConstraints", "Unit cell is infeasible." };
	}
}

inline void MathematicalCrystalChemistry::Design::CrystalDesigner::forceUpdateConstraints(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	constrainingCrystalStructure.reduceStructure();
	{
		if (constrainingCrystalStructure.hasFeasibleUnitCell())
		{
			constrainingCrystalStructure.updateTracingIndexPairs();
			constrainingCrystalStructure.createInteratomicDistanceConstraints();
			constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();

			m_interatomicDistanceTrackerUsing = 0;
			m_unitCellUsing = 0;
		}

		else
			throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "forceUpdateConstraints", "Unit cell is infeasible despite the reduction." };
	}
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_DESIGN_CRYSTALDESIGNER_H
