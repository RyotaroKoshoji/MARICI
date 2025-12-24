#include "CrystallineConstraintManager.h"

#include <cmath>

#include "GeometricalConstraintParameters.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components::Internal;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

CrystallineConstraintManager::CrystallineConstraintManager() noexcept
	: CrystalStructure{}
	, _feasibleErrorRate{ 0.0 }
	, _exclusiveRadiusRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultMinimumExclusionDistanceRatio() }
	, _interatomicDistanceTracerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerCutoffRatio() }
	, _interatomicDistanceConstrainerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceConstrainerCutoffRatio() }
	, _constrainingIndexPairs{}
	, _tracingIndexPairs{}
{
}

CrystallineConstraintManager::CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell& cell) noexcept
	: CrystalStructure{ cell }
	, _feasibleErrorRate{ 0.0 }
	, _exclusiveRadiusRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultMinimumExclusionDistanceRatio() }
	, _interatomicDistanceTracerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerCutoffRatio() }
	, _interatomicDistanceConstrainerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceConstrainerCutoffRatio() }
	, _constrainingIndexPairs{}
	, _tracingIndexPairs{}
{
}

CrystallineConstraintManager::CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell& cell, const std::vector<ConstrainingAtom>& atoms) noexcept
	: CrystalStructure{ cell, atoms }
	, _feasibleErrorRate{ 0.0 }
	, _exclusiveRadiusRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultMinimumExclusionDistanceRatio() }
	, _interatomicDistanceTracerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerCutoffRatio() }
	, _interatomicDistanceConstrainerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceConstrainerCutoffRatio() }
	, _constrainingIndexPairs{}
	, _tracingIndexPairs{}
{
}

CrystallineConstraintManager::CrystallineConstraintManager(const ChemToolkit::Crystallography::UnitCell& cell, std::vector<ConstrainingAtom>&& atoms) noexcept
	: CrystalStructure{ cell, std::move(atoms) }
	, _feasibleErrorRate{ 0.0 }
	, _exclusiveRadiusRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultMinimumExclusionDistanceRatio() }
	, _interatomicDistanceTracerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceTracerCutoffRatio() }
	, _interatomicDistanceConstrainerCutoffRatio{ MathematicalCrystalChemistry::Design::Optimization::GeometricalConstraintParameters::defaultInteratomicDistanceConstrainerCutoffRatio() }
	, _constrainingIndexPairs{}
	, _tracingIndexPairs{}
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

void CrystallineConstraintManager::updateTracingIndexPairs()
{
	clearInteratomicDistanceConstraints();
	NumericalMatrix inverseBasisVectors = unitCell().getInverseBasisVectors();


	for (size_type originalIndex = 0; originalIndex < atoms().size(); ++originalIndex)
	{
		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (isIonicAttractive(originalIndex, translatedOriginalIndex))
			{
				double neighborZoneRadius = _interatomicDistanceTracerCutoffRatio * _exclusiveRadiusRatio * (atoms()[originalIndex].ionicRadius().maximum() + atoms()[translatedOriginalIndex].ionicRadius().maximum());
				{
					for (const auto& latticePoint : enumerateNeighborLatticePoints(originalIndex, translatedOriginalIndex, inverseBasisVectors, neighborZoneRadius))
					{
						if (!(isOriginalLatticePoint(latticePoint)))
						{
							if (isTraceableIonicExclusionDistance(originalIndex, translatedOriginalIndex, latticePoint))
								_tracingIndexPairs.push_back(ConstrainerIndices<TranslatedAtomIndex>{ OriginalAtomIndex{ originalIndex }, TranslatedAtomIndex{ translatedOriginalIndex, latticePoint } });
						}
					}
				}
			}

			else
			{
				if (isIonicRepulsive(originalIndex, translatedOriginalIndex))
				{
					double neighborZoneRadius = _interatomicDistanceTracerCutoffRatio * (atoms()[originalIndex].ionicRepulsionRadius().minimum() + atoms()[translatedOriginalIndex].ionicRepulsionRadius().minimum());
					{
						for (const auto& latticePoint : enumerateNeighborLatticePoints(originalIndex, translatedOriginalIndex, inverseBasisVectors, neighborZoneRadius))
						{
							if (!(isOriginalLatticePoint(latticePoint)))
							{
								if (isTraceableIonicRepulsionDistance(originalIndex, translatedOriginalIndex, latticePoint))
									_tracingIndexPairs.push_back(ConstrainerIndices<TranslatedAtomIndex>{ OriginalAtomIndex{ originalIndex }, TranslatedAtomIndex{ translatedOriginalIndex, latticePoint } });
							}
						}
					}
				}

				else
				{
					double neighborZoneRadius = _interatomicDistanceTracerCutoffRatio * _exclusiveRadiusRatio * (atoms()[originalIndex].covalentRadius().maximum() + atoms()[translatedOriginalIndex].covalentRadius().maximum());
					{
						for (const auto& latticePoint : enumerateNeighborLatticePoints(originalIndex, translatedOriginalIndex, inverseBasisVectors, neighborZoneRadius))
						{
							if (!(isOriginalLatticePoint(latticePoint)))
							{
								if (isTraceableCovalentExclusionDistance(originalIndex, translatedOriginalIndex, latticePoint))
									_tracingIndexPairs.push_back(ConstrainerIndices<TranslatedAtomIndex>{ OriginalAtomIndex{ originalIndex }, TranslatedAtomIndex{ translatedOriginalIndex, latticePoint } });
							}
						}
					}
				}
			}
		}
	}

	for (size_type originalIndex = 0; originalIndex < atoms().size(); ++originalIndex)
	{
		LatticePoint originalLatticePoint{ 0,0,0 };


		if (isIonicRepulsive(originalIndex, originalIndex))
		{
			double neighborZoneRadius = _interatomicDistanceTracerCutoffRatio * (atoms()[originalIndex].ionicRepulsionRadius().minimum() + atoms()[originalIndex].ionicRepulsionRadius().minimum());
			{
				for (const auto& latticePoint : enumerateNeighborLatticePoints(originalIndex, originalIndex, inverseBasisVectors, neighborZoneRadius))
				{
					if (originalLatticePoint < latticePoint)
					{
						if (isTraceableIonicRepulsionDistance(originalIndex, originalIndex, latticePoint))
							_tracingIndexPairs.push_back(ConstrainerIndices<TranslatedAtomIndex>{ OriginalAtomIndex{ originalIndex }, TranslatedAtomIndex{ originalIndex, latticePoint } });
					}
				}
			}
		}

		else
		{
			double neighborZoneRadius = _interatomicDistanceTracerCutoffRatio * _exclusiveRadiusRatio * (atoms()[originalIndex].covalentRadius().maximum() + atoms()[originalIndex].covalentRadius().maximum());
			{
				for (const auto& latticePoint : enumerateNeighborLatticePoints(originalIndex, originalIndex, inverseBasisVectors, neighborZoneRadius))
				{
					if (originalLatticePoint < latticePoint)
					{
						if (isTraceableCovalentExclusionDistance(originalIndex, originalIndex, latticePoint))
							_tracingIndexPairs.push_back(ConstrainerIndices<TranslatedAtomIndex>{ OriginalAtomIndex{ originalIndex }, TranslatedAtomIndex{ originalIndex, latticePoint } });
					}
				}
			}
		}
	}
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
// Protected methods

void CrystallineConstraintManager::updateConstrainingIndexPairs()
{
	_constrainingIndexPairs.clear();


	for (const auto& indices : _tracingIndexPairs)
	{
		if (isIonicAttractive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (isConstrainableIonicExclusionDistance(indices.originalAtomIndex(), indices.translatedAtomIndex()))
				_constrainingIndexPairs.push_back(indices);
		}

		else if (isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (isConstrainableIonicRepulsionDistance(indices.originalAtomIndex(), indices.translatedAtomIndex()))
				_constrainingIndexPairs.push_back(indices);
		}

		else
		{
			if (isConstrainableCovalentExclusionDistance(indices.originalAtomIndex(), indices.translatedAtomIndex()))
				_constrainingIndexPairs.push_back(indices);
		}
	}
}

// Protected methods
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

std::vector<CrystallineConstraintManager::LatticePoint> CrystallineConstraintManager::enumerateNeighborLatticePoints(const OriginalAtomIndex originalIndex, const OriginalAtomIndex translatedIndex, const NumericalMatrix& inverseBasisVectors, const double neighborZoneRadius) const
{
	std::vector<LatticePoint> latticePoints;
	{
		const NumericalVector fractionalCoordinate = inverseBasisVectors * atoms()[originalIndex].cartesianCoordinate();

		double lengthA = neighborZoneRadius;
		{
			double normSquare = inverseBasisVectors(0, 0) * inverseBasisVectors(0, 0);
			normSquare += inverseBasisVectors(0, 1) * inverseBasisVectors(0, 1);
			normSquare += inverseBasisVectors(0, 2) * inverseBasisVectors(0, 2);

			lengthA *= std::sqrt(normSquare);
		}

		double lengthB = neighborZoneRadius;
		{
			double normSquare = inverseBasisVectors(1, 0) * inverseBasisVectors(1, 0);
			normSquare += inverseBasisVectors(1, 1) * inverseBasisVectors(1, 1);
			normSquare += inverseBasisVectors(1, 2) * inverseBasisVectors(1, 2);

			lengthB *= std::sqrt(normSquare);
		}

		double lengthC = neighborZoneRadius;
		{
			double normSquare = inverseBasisVectors(2, 0) * inverseBasisVectors(2, 0);
			normSquare += inverseBasisVectors(2, 1) * inverseBasisVectors(2, 1);
			normSquare += inverseBasisVectors(2, 2) * inverseBasisVectors(2, 2);

			lengthC *= std::sqrt(normSquare);
		}


		short minA = static_cast<short>(std::floor(fractionalCoordinate[0] - lengthA));
		short maxA = static_cast<short>(std::floor(fractionalCoordinate[0] + lengthA));

		short minB = static_cast<short>(std::floor(fractionalCoordinate[1] - lengthB));
		short maxB = static_cast<short>(std::floor(fractionalCoordinate[1] + lengthB));

		short minC = static_cast<short>(std::floor(fractionalCoordinate[2] - lengthC));
		short maxC = static_cast<short>(std::floor(fractionalCoordinate[2] + lengthC));


		for (short indexA = minA; indexA <= maxA; ++indexA)
		{
			for (short indexB = minB; indexB <= maxB; ++indexB)
			{
				for (short indexC = minC; indexC <= maxC; ++indexC)
					latticePoints.push_back(LatticePoint{ indexA, indexB, indexC });
			}
		}
	}

	return latticePoints;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
