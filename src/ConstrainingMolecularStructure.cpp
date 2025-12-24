#include "ConstrainingMolecularStructure.h"

#include <cmath>
#include <algorithm>
#include <random>

#include "SpaceGroupTypeSearcher.h"

#include "ObjectiveMolecularStructure.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

ConstrainingMolecularStructure::ConstrainingMolecularStructure() noexcept
	: LinkedPolyhedraRetriever{}
	, m_randomEngine{}
{
	std::random_device seedGen;
	m_randomEngine = std::mt19937{ seedGen() };
}

ConstrainingMolecularStructure::ConstrainingMolecularStructure(const ObjectiveMolecularStructure& structure)
	: LinkedPolyhedraRetriever{ structure.unitCell() }
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
	}
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

void ConstrainingMolecularStructure::initialize() noexcept
{
	CrystalStructure::initialize();
}

void ConstrainingMolecularStructure::initialize(const std::vector<ConstrainingAtom>& atomicArrangement) noexcept
{
	atoms() = atomicArrangement;
	centerAtoms();
	CrystallineConstraintManager::clearInteratomicDistanceConstraints();

	unitCell().basisVectors() *= std::pow((20.0 * getAtomicSphereVolume()), 0.333333333333333);
}

void ConstrainingMolecularStructure::initialize(std::vector<ConstrainingAtom>&& atomicArrangement) noexcept
{
	atoms() = std::move(atomicArrangement);
	centerAtoms();
	CrystallineConstraintManager::clearInteratomicDistanceConstraints();

	unitCell().basisVectors() *= std::pow((20.0 * getAtomicSphereVolume()), 0.333333333333333);
}

void ConstrainingMolecularStructure::importStructure(const ObjectiveMolecularStructure& structure)
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

void ConstrainingMolecularStructure::import(const ObjectiveMolecularStructure& structure)
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
	}
}

void ConstrainingMolecularStructure::distortStructure()
{
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
	}
}

void ConstrainingMolecularStructure::distortStructureLargely()
{
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
	}
}

void ConstrainingMolecularStructure::eraseInfeasibleIonicPolyhedraConnections()
{
	for (const auto& ionicCoordinationPolyhedraLinking : getCoordinationPolyhedraLinkings())
	{
		std::pair<IonicAtomicNumber, IonicAtomicNumber> ionicAtomicNumberKey = toIonicAtomicNumberPair(ionicCoordinationPolyhedraLinking.first);
		ChemicalComposition commonBridgingAnions = toChemicalComposition(ionicCoordinationPolyhedraLinking.second);
		ChemicalComposition closestFeasibleCommonBridgingAnions = feasiblePolyhedraConnections().getClosestFeasibleCommonBridgingComposition(ionicAtomicNumberKey.first, ionicAtomicNumberKey.second, commonBridgingAnions);


		if (closestFeasibleCommonBridgingAnions < commonBridgingAnions)
		{
			for (const auto& numAndCount : closestFeasibleCommonBridgingAnions)
			{
				std::vector<TranslatedAtomIndex> commonBridgingIndices;
				{
					for (const auto& translatedIndex : ionicCoordinationPolyhedraLinking.second)
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
						if (willChooseOriginalAtomIndex(ionicCoordinationPolyhedraLinking.first))
							eraseIonicBond(ionicCoordinationPolyhedraLinking.first.originalAtomIndex(), commonBridgingIndices.back().originalIndex());

						else
						{
							if (ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().isInOriginalCell())
								eraseIonicBond(ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back().originalIndex());

							else
							{
								TranslatedAtomIndex reverseIndex{ commonBridgingIndices.back().originalIndex(), ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().latticePoint() };
								reverseIndex.reverseLatticePoint();
								eraseIonicBond(ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), reverseIndex);
							}
						}
					}

					else
					{
						if (willChooseOriginalAtomIndex(ionicCoordinationPolyhedraLinking.first))
							eraseIonicBond(ionicCoordinationPolyhedraLinking.first.originalAtomIndex(), commonBridgingIndices.back());

						else
						{
							if (ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().isInOriginalCell())
								eraseIonicBond(ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back());

							else
							{
								TranslatedAtomIndex reverseIndex{ commonBridgingIndices.back().originalIndex(), ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().getRelativeLatticePointTo(commonBridgingIndices.back()) };

								if (reverseIndex.isInOriginalCell())
									eraseIonicBond(ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), commonBridgingIndices.back().originalIndex());
								else
									eraseIonicBond(ionicCoordinationPolyhedraLinking.first.translatedAtomIndex().originalIndex(), reverseIndex);
							}
						}
					}

					commonBridgingIndices.pop_back();
				}
			}
		}
	}
}

bool ConstrainingMolecularStructure::isFeasible() const
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

void ConstrainingMolecularStructure::createChemicalBonds()
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
}

void ConstrainingMolecularStructure::optimizeCoordinationCompositions()
{
	for (size_type centralAtomIndex = 0; centralAtomIndex < atoms().size(); ++centralAtomIndex)
	{
		const ConstrainingAtom& centralAtom = atoms()[centralAtomIndex];


		if (centralAtom.coordinationConstraints().hasFeasibleCompositions())
			makeClosestCoordinationComposition(centralAtomIndex);
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

void ConstrainingMolecularStructure::makeClosestCoordinationComposition(const OriginalAtomIndex centralAtomIndex)
{
	const ConstrainingAtom& centralAtom = atoms()[centralAtomIndex];


	if (centralAtom.coordinationConstraints().maxCoordinationNumber() < centralAtom.getCoordinationNumber())
	{
		for (const auto& numAndCount : centralAtom.coordinationConstraints().getClosestFeasibleComposition(getCoordinationComposition(centralAtomIndex)))
		{
			std::vector<std::pair<double, OriginalAtomIndex>> originalAtomIndices = getOrderedChemicalBondedOriginalIndices(centralAtomIndex, numAndCount.first);


			while (numAndCount.second < (originalAtomIndices.size()))
			{
				eraseCovalentBond(centralAtomIndex, originalAtomIndices.back().second);
				eraseIonicBond(centralAtomIndex, originalAtomIndices.back().second);
				originalAtomIndices.pop_back();
			}
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
