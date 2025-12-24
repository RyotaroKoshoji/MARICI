#include "CrystalDesigner.h"

using namespace MathematicalCrystalChemistry::Design;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

CrystalDesigner::CrystalDesigner() noexcept
	: _maxTotalStructuralOptimizing{ 0 }
	, _maxCeaselessGlobalStructuralOptimizing{ 0 }
	, _geometricalConstraintParameters{}
	, _randomStructureGenerator{}
	, _globalStructuralOptimizer{}
	, _localStructuralOptimizer{}
	, _preciseStructuralOptimizer{}
	, m_totalStructuralOptimizing{ 0 }
	, m_ceaselessGlobalStructuralOptimizing{ 0 }
	, m_interatomicDistanceTrackerUsing{ 0 }
	, m_unitCellUsing{ 0 }
{
}

CrystalDesigner::CrystalDesigner(const CrystalDesignParameters& parameters)
	: _maxTotalStructuralOptimizing{ parameters.maxTotalStructuralOptimizing() }
	, _maxCeaselessGlobalStructuralOptimizing{ parameters.maxCeaselessGlobalStructuralOptimizing() }
	, _geometricalConstraintParameters{ parameters.geometricalConstraintParameters() }
	, _randomStructureGenerator{ parameters.initialStructureGenerationParameters().randomStructureGenerationParameters() }
	, _globalStructuralOptimizer{ parameters.globalStructuralOptimizationParameters(), parameters.geometricalConstraintParameters() }
	, _localStructuralOptimizer{ parameters.localStructuralOptimizationParameters(), parameters.geometricalConstraintParameters() }
	, _preciseStructuralOptimizer{ parameters.preciseStructuralOptimizationParameters(), parameters.geometricalConstraintParameters() }
	, m_totalStructuralOptimizing{ 0 }
	, m_ceaselessGlobalStructuralOptimizing{ 0 }
	, m_interatomicDistanceTrackerUsing{ 0 }
	, m_unitCellUsing{ 0 }
{
}

// Constructors
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

void CrystalDesigner::execute(ConstrainingCrystalStructure& constrainingCrystalStructure) const
{
	initializeTimers();
	constrainingCrystalStructure.setFeasibleErrorRate(_globalStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());
	constrainingCrystalStructure.setExclusiveRadiusRatio(_geometricalConstraintParameters.minimumExclusionDistanceRatio());
	constrainingCrystalStructure.setInteratomicDistanceTracerCutoffRatio(_geometricalConstraintParameters.interatomicDistanceTracerCutoffRatio());
	constrainingCrystalStructure.setInteratomicDistanceConstrainerCutoffRatio(_geometricalConstraintParameters.interatomicDistanceConstrainerCutoffRatio());
	constrainingCrystalStructure.updateTracingIndexPairs();
	constrainingCrystalStructure.createInteratomicDistanceConstraints();
	constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();


	while (m_totalStructuralOptimizing < _maxTotalStructuralOptimizing)
	{
		constrainingCrystalStructure.setFeasibleErrorRate(_globalStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());

		ObjectiveCrystalStructure objectiveCrystalStructure{ constrainingCrystalStructure };
		applyGlobalStructuralOptimization(objectiveCrystalStructure);

		constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
		updateConstraints(constrainingCrystalStructure);


		if (constrainingCrystalStructure.isFeasibleCoordinationComposition())
		{
			m_ceaselessGlobalStructuralOptimizing = 0;
			objectiveCrystalStructure.import(constrainingCrystalStructure);


			if (applyLocalStructuralOptimization(objectiveCrystalStructure))
			{
				if (applyPreciseStructuralOptimization(objectiveCrystalStructure))
				{
					constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
					constrainingCrystalStructure.setFeasibleErrorRate(_preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());


					if (isFeasible(constrainingCrystalStructure))
					{
						initializeTimers();
						return;
					}

					else
					{
						constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
						constrainingCrystalStructure.distortStructure();
						m_unitCellUsing = _geometricalConstraintParameters.unitCellReductionTimeout();

						continue;
					}
				}
			}

			else
			{
				constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
				constrainingCrystalStructure.setFeasibleErrorRate(_localStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());

				reduceStructure(constrainingCrystalStructure);

				constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
				constrainingCrystalStructure.distortStructure();
				m_unitCellUsing = _geometricalConstraintParameters.unitCellReductionTimeout();

				continue;
			}
		}

		else
		{
			if (_maxCeaselessGlobalStructuralOptimizing < m_ceaselessGlobalStructuralOptimizing)
			{
				constrainingCrystalStructure.distortStructureLargely();
				reduceStructure(constrainingCrystalStructure);

				m_ceaselessGlobalStructuralOptimizing = 0;
			}
		}
	}


	constrainingCrystalStructure.setFeasibleErrorRate(_preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());
	constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
	initializeTimers();
}

void CrystalDesigner::execute(ConstrainingCrystalStructure& constrainingCrystalStructure, CrystalDesignRecorder& crystalDesignRecorder) const
{
	initializeTimers();
	constrainingCrystalStructure.setFeasibleErrorRate(_globalStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());
	constrainingCrystalStructure.setExclusiveRadiusRatio(_geometricalConstraintParameters.minimumExclusionDistanceRatio());
	constrainingCrystalStructure.setInteratomicDistanceTracerCutoffRatio(_geometricalConstraintParameters.interatomicDistanceTracerCutoffRatio());
	constrainingCrystalStructure.setInteratomicDistanceConstrainerCutoffRatio(_geometricalConstraintParameters.interatomicDistanceConstrainerCutoffRatio());
	constrainingCrystalStructure.updateTracingIndexPairs();
	constrainingCrystalStructure.createInteratomicDistanceConstraints();
	constrainingCrystalStructure.eraseInfeasibleIonicPolyhedraConnections();

	crystalDesignRecorder.forceRecord(ObjectiveCrystalStructure{ constrainingCrystalStructure });


	while (m_totalStructuralOptimizing < _maxTotalStructuralOptimizing)
	{
		constrainingCrystalStructure.setFeasibleErrorRate(_globalStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());
		constrainingCrystalStructure.normalizeAverageFractionalCoordinates();

		ObjectiveCrystalStructure objectiveCrystalStructure{ constrainingCrystalStructure };
		applyGlobalStructuralOptimization(objectiveCrystalStructure, crystalDesignRecorder);

		constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
		updateConstraints(constrainingCrystalStructure);


		if (constrainingCrystalStructure.isFeasibleCoordinationComposition())
		{
			m_ceaselessGlobalStructuralOptimizing = 0;
			objectiveCrystalStructure.import(constrainingCrystalStructure);


			if (applyLocalStructuralOptimization(objectiveCrystalStructure, crystalDesignRecorder))
			{
				if (applyPreciseStructuralOptimization(objectiveCrystalStructure, crystalDesignRecorder))
				{
					constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
					constrainingCrystalStructure.setFeasibleErrorRate(_preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());


					if (isFeasible(constrainingCrystalStructure))
					{
						initializeTimers();
						return;
					}

					else
					{
						constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
						constrainingCrystalStructure.distortStructure();
						m_unitCellUsing = _geometricalConstraintParameters.unitCellReductionTimeout();

						continue;
					}
				}
			}

			else
			{
				constrainingCrystalStructure.importStructure(objectiveCrystalStructure);
				constrainingCrystalStructure.setFeasibleErrorRate(_localStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());

				reduceStructure(constrainingCrystalStructure);

				constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
				constrainingCrystalStructure.distortStructure();
				m_unitCellUsing = _geometricalConstraintParameters.unitCellReductionTimeout();

				continue;
			}
		}

		else
		{
			if (_maxCeaselessGlobalStructuralOptimizing < m_ceaselessGlobalStructuralOptimizing)
			{
				constrainingCrystalStructure.distortStructureLargely();
				reduceStructure(constrainingCrystalStructure);

				m_ceaselessGlobalStructuralOptimizing = 0;
			}
		}
	}

		
	constrainingCrystalStructure.setFeasibleErrorRate(_preciseStructuralOptimizer.structuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());
	constrainingCrystalStructure.eraseInfeasibleChemicalBonds();
	initializeTimers();
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
