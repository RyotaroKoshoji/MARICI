#include "CoordinationPolyhedraRetriever.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components::Internal;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

CoordinationPolyhedraRetriever::CoordinationPolyhedraRetriever() noexcept
	: CrystallineConstraintManager{}
{
}

CoordinationPolyhedraRetriever::CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell& cell) noexcept
	: CrystallineConstraintManager{ cell }
{
}

CoordinationPolyhedraRetriever::CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell& cell, const std::vector<ConstrainingAtom>& atoms) noexcept
	: CrystallineConstraintManager{ cell, atoms }
{
}

CoordinationPolyhedraRetriever::CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell& cell, std::vector<ConstrainingAtom>&& atoms) noexcept
	: CrystallineConstraintManager{ cell, std::move(atoms) }
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
// Public methods

void CoordinationPolyhedraRetriever::eraseInfeasibleChemicalBonds()
{
	for (size_type originalIndex = 0; originalIndex < atoms().size(); ++originalIndex)
	{
		const ConstrainingAtom& constrainingAtom = atoms()[originalIndex];


		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (constrainingAtom.hasIonicBondWith(translatedOriginalIndex))
			{
				if (!(isFeasibleIonicBond(originalIndex, translatedOriginalIndex)))
					eraseIonicBond(originalIndex, translatedOriginalIndex);
			}

			if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
			{
				if (!(isFeasibleCovalentBond(originalIndex, translatedOriginalIndex)))
					eraseCovalentBond(originalIndex, translatedOriginalIndex);
			}
		}
	}

	for (const auto& indices : constrainingIndexPairs())
	{
		const ConstrainingAtom& constrainingAtom = atoms()[indices.originalAtomIndex()];


		if (constrainingAtom.hasIonicBondWith(indices.translatedAtomIndex()))
		{
			if (!(isFeasibleIonicBond(indices.originalAtomIndex(), indices.translatedAtomIndex())))
				eraseIonicBond(indices.originalAtomIndex(), indices.translatedAtomIndex());
		}

		if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
		{
			if (!(isFeasibleCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex())))
				eraseCovalentBond(indices.originalAtomIndex(), indices.translatedAtomIndex());
		}
	}
}

// Public methods
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

bool CoordinationPolyhedraRetriever::hasFeasibleCoordinationComposition(const OriginalAtomIndex originalAtomIndex) const
{
	const ConstrainingAtom& atom = atoms()[originalAtomIndex];


	if (atom.coordinationConstraints().hasFeasibleCompositions())
		return !(atom.coordinationConstraints().isInfeasibleComposition(getCoordinationComposition(originalAtomIndex)));

	else
	{
		if (atom.coordinationConstraints().hasFeasibleCovalentCoordinationNumbers())
		{
			if (atom.coordinationConstraints().isInfeasibleCovalentCoordinationNumber(atom.getCovalentCoordinationNumber()))
				return false;
		}

		if (atom.coordinationConstraints().hasFeasibleIonicCoordinationNumbers())
		{
			if (atom.coordinationConstraints().isInfeasibleIonicCoordinationNumber(atom.getIonicCoordinationNumber()))
				return false;
		}

		if (atom.coordinationConstraints().hasLowerBoundCompositions())
		{
			if (!(atom.coordinationConstraints().satisfyLowerBound(getCoordinationComposition(originalAtomIndex))))
				return false;
		}


		return true;
	}
}

bool CoordinationPolyhedraRetriever::hasInfeasibleChemicalBonds(const OriginalAtomIndex originalAtomIndex, const double errorRate) const
{
	const ConstrainingAtom& centralAtom = atoms()[originalAtomIndex];


	for (const auto& originalCoordinatedIndex : centralAtom.getIonicBondedOriginalAtomIndices())
	{
		if (!(isFeasibleIonicBond(originalAtomIndex, originalCoordinatedIndex, errorRate)))
			return true;
	}

	for (const auto& originalCoordinatedIndex : centralAtom.getCovalentBondedOriginalAtomIndices())
	{
		if (!(isFeasibleCovalentBond(originalAtomIndex, originalCoordinatedIndex, errorRate)))
			return true;
	}

	for (const auto& translatedCoordinatedIndex : centralAtom.getIonicBondedTranslatedAtomIndices())
	{
		if (!(isFeasibleIonicBond(originalAtomIndex, translatedCoordinatedIndex, errorRate)))
			return true;
	}

	for (const auto& translatedCoordinatedIndex : centralAtom.getCovalentBondedTranslatedAtomIndices())
	{
		if (!(isFeasibleCovalentBond(originalAtomIndex, translatedCoordinatedIndex, errorRate)))
			return true;
	}


	return false;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::OriginalAtomIndex>> CoordinationPolyhedraRetriever::getOrderedCovalentBondedOriginalIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, OriginalAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto coordinatedIndex : centralAtom.getCovalentBondedOriginalAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex].cartesianCoordinate() - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, OriginalAtomIndex>& a, const std::pair<double, OriginalAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::TranslatedAtomIndex>> CoordinationPolyhedraRetriever::getOrderedCovalentBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, TranslatedAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto& coordinatedIndex : centralAtom.getCovalentBondedTranslatedAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, TranslatedAtomIndex>& a, const std::pair<double, TranslatedAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::OriginalAtomIndex>> CoordinationPolyhedraRetriever::getOrderedIonicBondedOriginalIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, OriginalAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto coordinatedIndex : centralAtom.getIonicBondedOriginalAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex].cartesianCoordinate() - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, OriginalAtomIndex>& a, const std::pair<double, OriginalAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::TranslatedAtomIndex>> CoordinationPolyhedraRetriever::getOrderedIonicBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, TranslatedAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto& coordinatedIndex : centralAtom.getIonicBondedTranslatedAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, TranslatedAtomIndex>& a, const std::pair<double, TranslatedAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::OriginalAtomIndex>> CoordinationPolyhedraRetriever::getOrderedChemicalBondedOriginalIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, OriginalAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto coordinatedIndex : centralAtom.getCovalentBondedOriginalAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex].cartesianCoordinate() - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		for (const auto coordinatedIndex : centralAtom.getIonicBondedOriginalAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex].cartesianCoordinate() - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, OriginalAtomIndex>& a, const std::pair<double, OriginalAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::TranslatedAtomIndex>> CoordinationPolyhedraRetriever::getOrderedChemicalBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const
{
	std::vector<std::pair<double, TranslatedAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto& coordinatedIndex : centralAtom.getCovalentBondedTranslatedAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		for (const auto& coordinatedIndex : centralAtom.getIonicBondedTranslatedAtomIndices())
		{
			NumericalVector displacement = atoms()[coordinatedIndex.originalIndex()].cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
			indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, TranslatedAtomIndex>& a, const std::pair<double, TranslatedAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::OriginalAtomIndex>> CoordinationPolyhedraRetriever::getOrderedChemicalBondedOriginalIndices(const OriginalAtomIndex centralIndex, const AtomicNumber coordinatedAtomicNumber) const
{
	std::vector<std::pair<double, OriginalAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto coordinatedIndex : centralAtom.getCovalentBondedOriginalAtomIndices())
		{
			const ConstrainingAtom& coordinatedAtom = atoms()[coordinatedIndex];

			if (coordinatedAtomicNumber == coordinatedAtom.ionicAtomicNumber().atomicNumber())
			{
				NumericalVector displacement = coordinatedAtom.cartesianCoordinate() - centralAtom.cartesianCoordinate();
				indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
			}
		}

		for (const auto coordinatedIndex : centralAtom.getIonicBondedOriginalAtomIndices())
		{
			const ConstrainingAtom& coordinatedAtom = atoms()[coordinatedIndex];

			if (coordinatedAtomicNumber == coordinatedAtom.ionicAtomicNumber().atomicNumber())
			{
				NumericalVector displacement = coordinatedAtom.cartesianCoordinate() - centralAtom.cartesianCoordinate();
				indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
			}
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, OriginalAtomIndex>& a, const std::pair<double, OriginalAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<std::pair<double, CoordinationPolyhedraRetriever::TranslatedAtomIndex>> CoordinationPolyhedraRetriever::getOrderedChemicalBondedTranslatedIndices(const OriginalAtomIndex centralIndex, const AtomicNumber coordinatedAtomicNumber) const
{
	std::vector<std::pair<double, TranslatedAtomIndex>> indices;
	{
		const ConstrainingAtom& centralAtom = atoms()[centralIndex];

		for (const auto& coordinatedIndex : centralAtom.getCovalentBondedTranslatedAtomIndices())
		{
			const ConstrainingAtom& coordinatedAtom = atoms()[coordinatedIndex.originalIndex()];

			if (coordinatedAtomicNumber == coordinatedAtom.ionicAtomicNumber().atomicNumber())
			{
				NumericalVector displacement = coordinatedAtom.cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
				indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
			}
		}

		for (const auto& coordinatedIndex : centralAtom.getIonicBondedTranslatedAtomIndices())
		{
			const ConstrainingAtom& coordinatedAtom = atoms()[coordinatedIndex.originalIndex()];

			if (coordinatedAtomicNumber == coordinatedAtom.ionicAtomicNumber().atomicNumber())
			{
				NumericalVector displacement = coordinatedAtom.cartesianCoordinate() + toTranslationVector(coordinatedIndex.latticePoint()) - centralAtom.cartesianCoordinate();
				indices.push_back(std::make_pair(displacement.normSquare(), coordinatedIndex));
			}
		}

		std::sort(indices.begin(), indices.end(), [](const std::pair<double, TranslatedAtomIndex>& a, const std::pair<double, TranslatedAtomIndex>& b) { return (a.first < b.first); });
	}

	return indices;
}

std::vector<CoordinationPolyhedraRetriever::OriginalAtomIndex> CoordinationPolyhedraRetriever::getCovalentExcludedOriginalIndices(const OriginalAtomIndex originalAtomIndex) const
{
	std::vector<OriginalAtomIndex> covalentExcludedOriginalAtomIndices;
	{
		for (size_type translatedOriginalIndex = 0; translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (!(isIonicAttractive(originalAtomIndex, translatedOriginalIndex)) && !(isIonicRepulsive(originalAtomIndex, translatedOriginalIndex)))
			{
				if (!(atoms()[originalAtomIndex].hasCovalentBondWith(translatedOriginalIndex)))
				{
					if (isConstrainableCovalentExclusionDistance(originalAtomIndex, translatedOriginalIndex))
						covalentExcludedOriginalAtomIndices.push_back(translatedOriginalIndex);
				}
			}
		}
	}

	return covalentExcludedOriginalAtomIndices;
}

std::vector<CoordinationPolyhedraRetriever::TranslatedAtomIndex> CoordinationPolyhedraRetriever::getCovalentExcludedTranslatedIndices(const OriginalAtomIndex originalAtomIndex) const
{
	std::vector<TranslatedAtomIndex> covalentExcludedTranslatedAtomIndices;
	{
		for (const auto& indexPair : constrainingIndexPairs())
		{
			if (indexPair.originalAtomIndex() == originalAtomIndex)
			{
				if (!(isIonicAttractive(originalAtomIndex, indexPair.translatedAtomIndex().originalIndex())) && !(isIonicRepulsive(originalAtomIndex, indexPair.translatedAtomIndex().originalIndex())))
				{
					if (!(atoms()[originalAtomIndex].hasCovalentBondWith(indexPair.translatedAtomIndex())))
						covalentExcludedTranslatedAtomIndices.push_back(indexPair.translatedAtomIndex());
				}
			}

			else
			{
				if (indexPair.translatedAtomIndex().originalIndex() == originalAtomIndex)
				{
					if (!(isIonicAttractive(indexPair.originalAtomIndex(), originalAtomIndex)) && !(isIonicRepulsive(indexPair.originalAtomIndex(), originalAtomIndex)))
					{
						TranslatedAtomIndex reverseIndex{ indexPair.originalAtomIndex(), indexPair.translatedAtomIndex().latticePoint() };
						reverseIndex.reverseLatticePoint();

						if (!(atoms()[originalAtomIndex].hasCovalentBondWith(reverseIndex)))
							covalentExcludedTranslatedAtomIndices.push_back(std::move(reverseIndex));
					}
				}
			}
		}
	}

	return covalentExcludedTranslatedAtomIndices;
}

std::vector<CoordinationPolyhedraRetriever::OriginalAtomIndex> CoordinationPolyhedraRetriever::getIonicExcludedOriginalIndices(const OriginalAtomIndex originalAtomIndex) const
{
	std::vector<OriginalAtomIndex> ionicExcludedOriginalAtomIndices;
	{
		for (size_type translatedOriginalIndex = 0; translatedOriginalIndex < atoms().size(); ++translatedOriginalIndex)
		{
			if (isIonicAttractive(originalAtomIndex, translatedOriginalIndex))
			{
				if (!(atoms()[originalAtomIndex].hasIonicBondWith(translatedOriginalIndex)))
				{
					if (isConstrainableIonicExclusionDistance(originalAtomIndex, translatedOriginalIndex))
						ionicExcludedOriginalAtomIndices.push_back(translatedOriginalIndex);
				}
			}
		}
	}

	return ionicExcludedOriginalAtomIndices;
}

std::vector<CoordinationPolyhedraRetriever::TranslatedAtomIndex> CoordinationPolyhedraRetriever::getIonicExcludedTranslatedIndices(const OriginalAtomIndex originalAtomIndex) const
{
	std::vector<TranslatedAtomIndex> ionicExcludedTranslatedAtomIndices;
	{
		for (const auto& indexPair : constrainingIndexPairs())
		{
			if (indexPair.originalAtomIndex() == originalAtomIndex)
			{
				if (isIonicAttractive(originalAtomIndex, indexPair.translatedAtomIndex().originalIndex()))
				{
					if (!(atoms()[originalAtomIndex].hasIonicBondWith(indexPair.translatedAtomIndex())))
						ionicExcludedTranslatedAtomIndices.push_back(indexPair.translatedAtomIndex());
				}
			}

			else
			{
				if (indexPair.translatedAtomIndex().originalIndex() == originalAtomIndex)
				{
					if (isIonicAttractive(indexPair.originalAtomIndex(), originalAtomIndex))
					{
						TranslatedAtomIndex reverseIndex{ indexPair.originalAtomIndex(), indexPair.translatedAtomIndex().latticePoint() };
						reverseIndex.reverseLatticePoint();

						if (!(atoms()[originalAtomIndex].hasIonicBondWith(reverseIndex)))
							ionicExcludedTranslatedAtomIndices.push_back(std::move(reverseIndex));
					}
				}
			}
		}
	}

	return ionicExcludedTranslatedAtomIndices;
}

// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
