#include "ConstrainingCrystalStructure.h"

#include <algorithm>
#include <random>

#include "SpaceGroupTypeSearcher.h"

#include "ObjectiveCrystalStructure.h"
#include "OptimalCrystalStructure.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

double ConstrainingCrystalStructure::s_defaultMinimumFeasiblePackingFraction{ 0.2 };



ConstrainingCrystalStructure::ConstrainingCrystalStructure() noexcept
	: LinkedPolyhedraRetriever{}
	, _minimumFeasiblePackingFraction{ s_defaultMinimumFeasiblePackingFraction }
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };
}

ConstrainingCrystalStructure::ConstrainingCrystalStructure(const ChemToolkit::Crystallography::UnitCell& cell, const std::vector<ConstrainingAtom>& atoms) noexcept
	: LinkedPolyhedraRetriever{ cell, atoms }
	, _minimumFeasiblePackingFraction{ s_defaultMinimumFeasiblePackingFraction }
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };
}

ConstrainingCrystalStructure::ConstrainingCrystalStructure(const ChemToolkit::Crystallography::UnitCell& cell, std::vector<ConstrainingAtom>&& atoms) noexcept
	: LinkedPolyhedraRetriever{ cell, std::move(atoms) }
	, _minimumFeasiblePackingFraction{ s_defaultMinimumFeasiblePackingFraction }
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };
}

ConstrainingCrystalStructure::ConstrainingCrystalStructure(const ObjectiveCrystalStructure& structure)
	: LinkedPolyhedraRetriever{ structure.unitCell() }
	, _minimumFeasiblePackingFraction{ s_defaultMinimumFeasiblePackingFraction }
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };

	{
		if (structure.isValid())
		{
			for (size_type index = 0; index < structure.atoms().size(); ++index)
				atoms().push_back(ConstrainingAtom{ structure.correspondingIonicAtomicNumbers()[index], structure.correspondingCoordinationConstraints()[index], structure.atoms()[index] });
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "import", "Argument objective crystal structure is invalid." };
	}


	updateTracingIndexPairs();
	updateConstrainingIndexPairs();
	{
		for (const auto& originalConstrainerIndices : structure.covalentBondedIndices())
			createCovalentBond(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& originalConstrainerIndices : structure.ionicBondedIndices())
			createIonicBond(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& originalConstrainerIndices : structure.ionicRepulsedIndices())
			createIonicRepulsion(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedCovalentBondedIndices())
			createCovalentBond(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedIonicBondedIndices())
			createIonicBond(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedIonicRepulsedIndices())
			createIonicRepulsion(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());
	}
}

ConstrainingCrystalStructure::ConstrainingCrystalStructure(const OptimalCrystalStructure& structure)
	: LinkedPolyhedraRetriever{ structure.unitCell() }
	, _minimumFeasiblePackingFraction{ s_defaultMinimumFeasiblePackingFraction }
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };


	for (const auto& atom : structure.atoms())
		atoms().push_back(ConstrainingAtom{ atom });
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

void ConstrainingCrystalStructure::importStructure(const ObjectiveCrystalStructure& structure)
{
	unitCell() = structure.unitCell();
	atoms().clear();
	{
		if (structure.isValid())
		{
			for (size_type index = 0; index < structure.atoms().size(); ++index)
				atoms().push_back(ConstrainingAtom{ structure.correspondingIonicAtomicNumbers()[index], structure.correspondingCoordinationConstraints()[index], structure.atoms()[index] });
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "import", "Argument objective crystal structure is invalid." };
	}
}

void ConstrainingCrystalStructure::import(const ObjectiveCrystalStructure& structure)
{
	unitCell() = structure.unitCell();
	atoms().clear();
	{
		if (structure.isValid())
		{
			for (size_type index = 0; index < structure.atoms().size(); ++index)
				atoms().push_back(ConstrainingAtom{ structure.correspondingIonicAtomicNumbers()[index], structure.correspondingCoordinationConstraints()[index], structure.atoms()[index] });
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "import", "Argument objective crystal structure is invalid." };
	}

	
	clearCovalentBonds();
	clearIonicBonds();
	clearIonicRepulsions();
	{
		for (const auto& originalConstrainerIndices : structure.covalentBondedIndices())
			createCovalentBond(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& originalConstrainerIndices : structure.ionicBondedIndices())
			createIonicBond(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& originalConstrainerIndices : structure.ionicRepulsedIndices())
			createIonicRepulsion(originalConstrainerIndices.originalAtomIndex(), originalConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedCovalentBondedIndices())
			createCovalentBond(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedIonicBondedIndices())
			createIonicBond(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());

		for (const auto& translatedConstrainerIndices : structure.translatedIonicRepulsedIndices())
			createIonicRepulsion(translatedConstrainerIndices.originalAtomIndex(), translatedConstrainerIndices.translatedAtomIndex());
	}
}

void ConstrainingCrystalStructure::import(const OptimalCrystalStructure& structure)
{
	unitCell().basisVectors() = structure.unitCell().basisVectors();
	atoms().clear();
	{
		for (const auto& atom : structure.atoms())
			atoms().push_back(ConstrainingAtom{ atom });
	}

	CrystallineConstraintManager::clearInteratomicDistanceConstraints();
}

void ConstrainingCrystalStructure::distortStructure()
{
	NumericalMatrix stressTensor;
	{
		std::uniform_real_distribution<double> stressTensorDistributor{ -0.1, 0.1 };

		stressTensor(0, 0) = stressTensorDistributor(m_randomEngine);
		stressTensor(1, 0) = stressTensorDistributor(m_randomEngine);
		stressTensor(2, 0) = stressTensorDistributor(m_randomEngine);

		stressTensor(0, 1) = stressTensor(1, 0);
		stressTensor(1, 1) = stressTensorDistributor(m_randomEngine);
		stressTensor(2, 1) = stressTensorDistributor(m_randomEngine);

		stressTensor(0, 2) = stressTensor(2, 0);
		stressTensor(1, 2) = stressTensor(2, 1);
		stressTensor(2, 2) = stressTensorDistributor(m_randomEngine);
	}

	{
		NumericalMatrix diplacementLatticeVectors = stressTensor * unitCell().basisVectors();
		unitCell().basisVectors() += diplacementLatticeVectors;
	}


	std::uniform_real_distribution<double> displacementDistributor{ 0.0, 0.1 };
	std::uniform_real_distribution<double> angleDistributor{ 0.0, MathToolkit::Generic::MathConstants::getPi() };

	for (auto& atom : atoms())
	{
		NumericalVector displacement;
		{
			double atomicRadius = 0.0;
			{
				if (atom.ionicAtomicNumber().formalCharge().isAnion() || atom.ionicAtomicNumber().formalCharge().isCation())
					atomicRadius = atom.ionicRadius().maximum();
				else
					atomicRadius = atom.covalentRadius().maximum();
			}

			double displacementSize = displacementDistributor(m_randomEngine) * atomicRadius;
			double theta = angleDistributor(m_randomEngine);
			double phi = 2.0 * angleDistributor(m_randomEngine);

			displacement[0] = displacementSize * std::sin(theta) * std::cos(phi);
			displacement[1] = displacementSize * std::sin(theta) * std::sin(phi);
			displacement[2] = displacementSize * std::cos(theta);
		}
		atom.cartesianCoordinate() += displacement;


		NumericalVector stressDisplacement = stressTensor * atom.cartesianCoordinate();
		atom.cartesianCoordinate() += stressDisplacement;
	}
}

void ConstrainingCrystalStructure::distortStructureLargely()
{
	NumericalMatrix stressTensor;
	{
		std::uniform_real_distribution<double> stressTensorDistributor{ -0.2, 0.2 };

		stressTensor(0, 0) = stressTensorDistributor(m_randomEngine);
		stressTensor(1, 0) = stressTensorDistributor(m_randomEngine);
		stressTensor(2, 0) = stressTensorDistributor(m_randomEngine);

		stressTensor(0, 1) = stressTensor(1, 0);
		stressTensor(1, 1) = stressTensorDistributor(m_randomEngine);
		stressTensor(2, 1) = stressTensorDistributor(m_randomEngine);

		stressTensor(0, 2) = stressTensor(2, 0);
		stressTensor(1, 2) = stressTensor(2, 1);
		stressTensor(2, 2) = stressTensorDistributor(m_randomEngine);
	}

	{
		NumericalMatrix diplacementLatticeVectors = stressTensor * unitCell().basisVectors();
		unitCell().basisVectors() += diplacementLatticeVectors;
	}


	std::uniform_real_distribution<double> displacementDistributor{ 0.0, 0.2 };
	std::uniform_real_distribution<double> angleDistributor{ 0.0, MathToolkit::Generic::MathConstants::getPi() };

	for (auto& atom : atoms())
	{
		NumericalVector displacement;
		{
			double atomicRadius = 0.0;
			{
				if (atom.ionicAtomicNumber().formalCharge().isAnion() || atom.ionicAtomicNumber().formalCharge().isCation())
					atomicRadius = atom.ionicRadius().maximum();
				else
					atomicRadius = atom.covalentRadius().maximum();
			}

			double displacementSize = displacementDistributor(m_randomEngine) * atomicRadius;
			double theta = angleDistributor(m_randomEngine);
			double phi = 2.0 * angleDistributor(m_randomEngine);

			displacement[0] = displacementSize * std::sin(theta) * std::cos(phi);
			displacement[1] = displacementSize * std::sin(theta) * std::sin(phi);
			displacement[2] = displacementSize * std::cos(theta);
		}
		atom.cartesianCoordinate() += displacement;


		NumericalVector stressDisplacement = stressTensor * atom.cartesianCoordinate();
		atom.cartesianCoordinate() += stressDisplacement;
	}
}

void ConstrainingCrystalStructure::reduceStructure()
{
	double cellVolume = unitCell().getCellVolume();
	{
		if (0.0 < cellVolume)
		{
			if ((_minimumFeasiblePackingFraction * cellVolume) < getAtomicSphereVolume())
			{
				ChemToolkit::Crystallography::Symmetry::SpaceGroupTypeSearcher spaceGroupTypeSearcher{ 0.000001 };
				spaceGroupTypeSearcher.applyDelaunayReduction(unitCell(), atoms());

				CrystallineConstraintManager::clearInteratomicDistanceConstraints();
			}

			else
				throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "reduceStructure", "Unit cell cannot be reduced due to small packing fraction." };
		}

		else
			throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "reduceStructure", "Unit cell cannot be reduced since volume of unit cell is not more than zero." };
	}
}

void ConstrainingCrystalStructure::eraseInfeasibleIonicPolyhedraConnections()
{
	for (const auto& coordinationPolyhedraLinking : getCoordinationPolyhedraLinkings())
	{
		std::pair<IonicAtomicNumber, IonicAtomicNumber> ionicAtomicNumberKey = toIonicAtomicNumberPair(coordinationPolyhedraLinking.first);
		ChemicalComposition commonBridgingAnions = toChemicalComposition(coordinationPolyhedraLinking.second);
		ChemicalComposition closestFeasibleCommonBridgingAnions = feasiblePolyhedraConnections().getClosestFeasibleCommonBridgingComposition(ionicAtomicNumberKey.first, ionicAtomicNumberKey.second, commonBridgingAnions);


		if (closestFeasibleCommonBridgingAnions < commonBridgingAnions)
		{
			for (const auto& numAndCount : closestFeasibleCommonBridgingAnions)
			{
				std::vector<TranslatedAtomIndex> commonBridgingIndices;
				{
					for (const auto& translatedIndex : coordinationPolyhedraLinking.second)
					{
						if (atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber() == numAndCount.first)
							commonBridgingIndices.push_back(translatedIndex);
					}

					std::shuffle(commonBridgingIndices.begin(), commonBridgingIndices.end(), m_randomEngine);
				}


				for (size_type rep = numAndCount.second; rep < commonBridgingAnions.count(numAndCount.first); ++rep)
				{
					if (commonBridgingIndices.back().isInOriginalCell())
					{
						if (willChooseOriginalAtomIndex(coordinationPolyhedraLinking.first))
							eraseIonicBond(coordinationPolyhedraLinking.first.originalAtomIndex(), commonBridgingIndices.back().originalIndex());

						else
						{
							if (coordinationPolyhedraLinking.first.translatedAtomIndex().isInOriginalCell())
								eraseIonicBond(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back().originalIndex());

							else
							{
								TranslatedAtomIndex reverseIndex{ commonBridgingIndices.back().originalIndex(), coordinationPolyhedraLinking.first.translatedAtomIndex().latticePoint() };
								reverseIndex.reverseLatticePoint();
								eraseIonicBond(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), reverseIndex);
							}
						}
					}

					else
					{
						if (willChooseOriginalAtomIndex(coordinationPolyhedraLinking.first))
							eraseIonicBond(coordinationPolyhedraLinking.first.originalAtomIndex(), commonBridgingIndices.back());

						else
						{
							if (coordinationPolyhedraLinking.first.translatedAtomIndex().isInOriginalCell())
								eraseIonicBond(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back());

							else
							{
								TranslatedAtomIndex reverseIndex{ commonBridgingIndices.back().originalIndex(), coordinationPolyhedraLinking.first.translatedAtomIndex().getRelativeLatticePointTo(commonBridgingIndices.back()) };

								if (reverseIndex.isInOriginalCell())
									eraseIonicBond(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back().originalIndex());
								else
									eraseIonicBond(coordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), reverseIndex);
							}
						}
					}

					commonBridgingIndices.pop_back();
				}
			}
		}
	}
}

bool ConstrainingCrystalStructure::isFeasible() const
{
	for (size_type index = 0; index < atoms().size(); ++index)
	{
		if (!(hasFeasibleCoordinationComposition(index)))
			return false;
	}


	for (size_type originalIndex = 0; originalIndex < atoms().size(); ++originalIndex)
	{
		const ConstrainingAtom& constrainingAtom = atoms()[originalIndex];


		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (isIonicAttractive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasIonicBondWith(translatedOriginalIndex))
				{
					if (!(isFeasibleIonicBond(originalIndex, translatedOriginalIndex)))
						return false;
				}

				else
				{
					if (!(isFeasibleIonicExclusion(originalIndex, translatedOriginalIndex)))
						return false;
				}
			}

			else if (isIonicRepulsive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
				{
					if (!(isFeasibleCovalentBond(originalIndex, translatedOriginalIndex)))
						return false;
				}

				else
				{
					if (!(isFeasibleIonicRepulsion(originalIndex, translatedOriginalIndex)))
						return false;
				}
			}

			else
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
				{
					if (!(isFeasibleCovalentBond(originalIndex, translatedOriginalIndex)))
						return false;
				}

				else
				{
					if (!(isFeasibleCovalentExclusion(originalIndex, translatedOriginalIndex)))
						return false;
				}
			}
		}
	}


	for (const auto& indices : constrainingIndexPairs())
	{
		const ConstrainingAtom& constrainingAtom = atoms()[indices.originalAtomIndex()];


		if (isIonicAttractive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasIonicBondWith(indices.translatedAtomIndex()))
			{
				if (!(isFeasibleIonicBond(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}

			else
			{
				if (!(isFeasibleIonicExclusion(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}
		}

		else if (isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
			{
				if (!(isFeasibleCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}

			else
			{
				if (!(isFeasibleIonicRepulsion(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}
		}

		else
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
			{
				if (!(isFeasibleCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}

			else
			{
				if (!(isFeasibleCovalentExclusion(indices.originalAtomIndex(), indices.translatedAtomIndex())))
					return false;
			}
		}
	}

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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private methods

void ConstrainingCrystalStructure::createChemicalBonds()
{
	for (size_type originalIndex = 0; originalIndex < atoms().size(); ++originalIndex)
	{
		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (isInnateChemicalBondable(originalIndex, translatedOriginalIndex))
			{
				if (isIonicAttractive(originalIndex, translatedOriginalIndex))
				{
					if (isInnateIonicBondable(originalIndex, translatedOriginalIndex) && isConstrainableIonicBondingDistance(originalIndex, translatedOriginalIndex))
						createIonicBond(originalIndex, translatedOriginalIndex);
				}

				else if (isIonicRepulsive(originalIndex, translatedOriginalIndex))
				{
					if (isInnateCovalentBondable(originalIndex, translatedOriginalIndex) && isConstrainableCovalentBondingDistance(originalIndex, translatedOriginalIndex))
						createCovalentBond(originalIndex, translatedOriginalIndex);
					else
						createIonicRepulsion(originalIndex, translatedOriginalIndex);
				}

				else
				{
					if (isInnateCovalentBondable(originalIndex, translatedOriginalIndex) && isConstrainableCovalentBondingDistance(originalIndex, translatedOriginalIndex))
						createCovalentBond(originalIndex, translatedOriginalIndex);
				}
			}

			else
			{
				if (isIonicRepulsive(originalIndex, translatedOriginalIndex))
					createIonicRepulsion(originalIndex, translatedOriginalIndex);
			}
		}
	}


	for (const auto& indices : constrainingIndexPairs())
	{
		NumericalVector translationVector = toTranslationVector(indices.translatedAtomIndex().latticePoint());


		if (isInnateChemicalBondable(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (isIonicAttractive(indices.originalAtomIndex(), indices.translatedAtomIndex()))
			{
				if (isInnateIonicBondable(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()) && isConstrainableIonicBondingDistance(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex(), translationVector))
					createIonicBond(indices.originalAtomIndex(), indices.translatedAtomIndex());
			}

			else if (isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex()))
			{
				if (isInnateCovalentBondable(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()) && isConstrainableCovalentBondingDistance(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex(), translationVector))
					createCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex());
				else
					createIonicRepulsion(indices.originalAtomIndex(), indices.translatedAtomIndex());
			}

			else
			{
				if (isInnateCovalentBondable(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()) && isConstrainableCovalentBondingDistance(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex(), translationVector))
					createCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex());
			}
		}

		else
		{
			if (isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex()))
				createIonicRepulsion(indices.originalAtomIndex(), indices.translatedAtomIndex());
		}
	}
}

void ConstrainingCrystalStructure::optimizeCoordinationCompositions()
{
	for (size_type centralAtomIndex = 0; centralAtomIndex < atoms().size(); ++centralAtomIndex)
	{
		const ConstrainingAtom& centralAtom = atoms()[centralAtomIndex];


		if (centralAtom.coordinationConstraints().hasFeasibleCompositions())
			makeClosestCoordinationComposition(centralAtomIndex);

		else
		{
			if (centralAtom.coordinationConstraints().hasFeasibleCovalentCoordinationNumbers())
			{
				std::vector<std::pair<double, OriginalAtomIndex>> covalentOriginalAtomIndices = getOrderedCovalentBondedOriginalIndices(centralAtomIndex);
				std::vector<std::pair<double, TranslatedAtomIndex>> covalentTranslatedAtomIndices = getOrderedCovalentBondedTranslatedIndices(centralAtomIndex);
				std::vector<std::pair<double, OriginalAtomIndex>> ionicOriginalAtomIndices = getOrderedIonicBondedOriginalIndices(centralAtomIndex);
				std::vector<std::pair<double, TranslatedAtomIndex>> ionicTranslatedAtomIndices = getOrderedIonicBondedTranslatedIndices(centralAtomIndex);


				if (centralAtom.coordinationConstraints().hasLowerBoundCompositions())
				{
					ChemicalComposition closestLowerBoundComposition = centralAtom.coordinationConstraints().getClosestLowerBoundComposition(toCoordinationComposition(covalentOriginalAtomIndices, covalentTranslatedAtomIndices), toCoordinationComposition(ionicOriginalAtomIndices, ionicTranslatedAtomIndices));


					if (centralAtom.coordinationConstraints().getClosestFeasibleCovalentCoordinationNumber(centralAtom.getCovalentCoordinationNumber()) < centralAtom.getCovalentCoordinationNumber())
					{
						size_type closestCoordinationNumber = centralAtom.coordinationConstraints().getClosestFeasibleCovalentCoordinationNumber(centralAtom.getCovalentCoordinationNumber());
						makeClosestCoordinationNumber(centralAtomIndex, closestCoordinationNumber, closestLowerBoundComposition, covalentOriginalAtomIndices, covalentTranslatedAtomIndices);
					}

					if (centralAtom.coordinationConstraints().getClosestFeasibleIonicCoordinationNumber(centralAtom.getIonicCoordinationNumber()) < centralAtom.getIonicCoordinationNumber())
					{
						size_type closestCoordinationNumber = centralAtom.coordinationConstraints().getClosestFeasibleIonicCoordinationNumber(centralAtom.getIonicCoordinationNumber());
						makeClosestCoordinationNumber(centralAtomIndex, closestCoordinationNumber, closestLowerBoundComposition, ionicOriginalAtomIndices, ionicTranslatedAtomIndices);
					}
				}

				else
				{
					if (centralAtom.coordinationConstraints().getClosestFeasibleCovalentCoordinationNumber(centralAtom.getCovalentCoordinationNumber()) < centralAtom.getCovalentCoordinationNumber())
					{
						size_type closestCoordinationNumber = centralAtom.coordinationConstraints().getClosestFeasibleCovalentCoordinationNumber(centralAtom.getCovalentCoordinationNumber());
						makeClosestCoordinationNumber(centralAtomIndex, closestCoordinationNumber, covalentOriginalAtomIndices, covalentTranslatedAtomIndices);
					}

					if (centralAtom.coordinationConstraints().getClosestFeasibleIonicCoordinationNumber(centralAtom.getIonicCoordinationNumber()) < centralAtom.getIonicCoordinationNumber())
					{
						size_type closestCoordinationNumber = centralAtom.coordinationConstraints().getClosestFeasibleIonicCoordinationNumber(centralAtom.getIonicCoordinationNumber());
						makeClosestCoordinationNumber(centralAtomIndex, closestCoordinationNumber, ionicOriginalAtomIndices, ionicTranslatedAtomIndices);
					}
				}
			}
		}
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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private utility

void ConstrainingCrystalStructure::makeClosestCoordinationComposition(const OriginalAtomIndex centralAtomIndex)
{
	const ConstrainingAtom& centralAtom = atoms()[centralAtomIndex];


	if (centralAtom.coordinationConstraints().maxCoordinationNumber() < centralAtom.getCoordinationNumber())
	{
		for (const auto& numAndCount : centralAtom.coordinationConstraints().getClosestFeasibleComposition(getCoordinationComposition(centralAtomIndex)))
		{
			std::vector<std::pair<double, OriginalAtomIndex>> originalAtomIndices = getOrderedChemicalBondedOriginalIndices(centralAtomIndex, numAndCount.first);
			std::vector<std::pair<double, TranslatedAtomIndex>> translatedAtomIndices = getOrderedChemicalBondedTranslatedIndices(centralAtomIndex, numAndCount.first);


			while (numAndCount.second < (originalAtomIndices.size() + translatedAtomIndices.size()))
			{
				double originalLongestDistanceSquare = 0.0;
				double translatedLongestDistanceSquare = 0.0;
				{
					if (!(originalAtomIndices.empty()))
						originalLongestDistanceSquare = originalAtomIndices.back().first;

					if(!(translatedAtomIndices.empty()))
						translatedLongestDistanceSquare = translatedAtomIndices.back().first;
				}


				if (translatedLongestDistanceSquare < originalLongestDistanceSquare)
				{
					eraseCovalentBond(centralAtomIndex, originalAtomIndices.back().second);
					eraseIonicBond(centralAtomIndex, originalAtomIndices.back().second);
					originalAtomIndices.pop_back();
				}

				else
				{
					eraseCovalentBond(centralAtomIndex, translatedAtomIndices.back().second);
					eraseIonicBond(centralAtomIndex, translatedAtomIndices.back().second);
					translatedAtomIndices.pop_back();
				}
			}
		}
	}
}

void ConstrainingCrystalStructure::makeClosestCoordinationNumber(const OriginalAtomIndex centralAtomIndex, const size_type closestCoordinationNumber, std::vector<std::pair<double, OriginalAtomIndex>>& originalAtomIndices, std::vector<std::pair<double, TranslatedAtomIndex>>& translatedAtomIndices)
{
	while (closestCoordinationNumber < (originalAtomIndices.size() + translatedAtomIndices.size()))
	{
		double originalLongestDistanceSquare = 0.0;
		double translatedLongestDistanceSquare = 0.0;
		{
			if (!(originalAtomIndices.empty()))
				originalLongestDistanceSquare = originalAtomIndices.back().first;

			if (!(translatedAtomIndices.empty()))
				translatedLongestDistanceSquare = translatedAtomIndices.back().first;
		}


		if (translatedLongestDistanceSquare < originalLongestDistanceSquare)
		{
			eraseCovalentBond(centralAtomIndex, originalAtomIndices.back().second);
			eraseIonicBond(centralAtomIndex, originalAtomIndices.back().second);
			originalAtomIndices.pop_back();
		}

		else
		{
			eraseCovalentBond(centralAtomIndex, translatedAtomIndices.back().second);
			eraseIonicBond(centralAtomIndex, translatedAtomIndices.back().second);
			translatedAtomIndices.pop_back();
		}
	}
}

void ConstrainingCrystalStructure::makeClosestCoordinationNumber(const OriginalAtomIndex centralAtomIndex, const size_type closestCoordinationNumber, const ChemicalComposition& closestLowerBoundComposition, std::vector<std::pair<double, OriginalAtomIndex>>& originalAtomIndices, std::vector<std::pair<double, TranslatedAtomIndex>>& translatedAtomIndices)
{
	const ConstrainingAtom& centralAtom = atoms()[centralAtomIndex];


	while (closestCoordinationNumber < (originalAtomIndices.size() + translatedAtomIndices.size()))
	{
		ChemicalComposition coordinationComposition = toCoordinationComposition(originalAtomIndices, translatedAtomIndices);


		if (translatedAtomIndices.back().first < originalAtomIndices.back().first)
		{
			AtomicNumber longestAtomicNumber = atoms()[originalAtomIndices.back().second].ionicAtomicNumber().atomicNumber();

			if (closestLowerBoundComposition.count(longestAtomicNumber) < coordinationComposition.count(longestAtomicNumber))
			{
				eraseCovalentBond(centralAtomIndex, originalAtomIndices.back().second);
				eraseIonicBond(centralAtomIndex, originalAtomIndices.back().second);
				originalAtomIndices.pop_back();
			}

			else
				std::rotate(originalAtomIndices.begin(), (originalAtomIndices.end() - 1), originalAtomIndices.end());
		}

		else
		{
			AtomicNumber longestAtomicNumber = atoms()[translatedAtomIndices.back().second.originalIndex()].ionicAtomicNumber().atomicNumber();

			if (closestLowerBoundComposition.count(longestAtomicNumber) < coordinationComposition.count(longestAtomicNumber))
			{
				eraseCovalentBond(centralAtomIndex, translatedAtomIndices.back().second);
				eraseIonicBond(centralAtomIndex, translatedAtomIndices.back().second);
				translatedAtomIndices.pop_back();
			}

			else
				std::rotate(translatedAtomIndices.begin(), (translatedAtomIndices.end() - 1), translatedAtomIndices.end());
		}
	}
}

// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
