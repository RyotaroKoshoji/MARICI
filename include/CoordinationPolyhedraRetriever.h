#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_COORDINATIONPOLYHEDRARETRIEVER_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_COORDINATIONPOLYHEDRARETRIEVER_H

#include <utility>

#include "AtomicNumber.h"
#include "ChemicalComposition.h"

#include "CrystallineConstraintManager.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			namespace Internal
			{
				class CoordinationPolyhedraRetriever :public CrystallineConstraintManager
				{
				protected:
					using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
					using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;

// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Constructors, destructor, and operators

				protected:
					CoordinationPolyhedraRetriever() noexcept;
					explicit CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&) noexcept;
					CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&, const std::vector<ConstrainingAtom>&) noexcept;
					CoordinationPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&, std::vector<ConstrainingAtom>&&) noexcept;

				public:
					virtual ~CoordinationPolyhedraRetriever() = default;

					CoordinationPolyhedraRetriever(const CoordinationPolyhedraRetriever&) = default;
					CoordinationPolyhedraRetriever(CoordinationPolyhedraRetriever&&) noexcept = default;
					CoordinationPolyhedraRetriever& operator=(const CoordinationPolyhedraRetriever&) = default;
					CoordinationPolyhedraRetriever& operator=(CoordinationPolyhedraRetriever&&) noexcept = default;

				// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Public methods

					void eraseInfeasibleChemicalBonds();

				// Public methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Protected methods

				protected:
					void clearChemicalBonds(const OriginalAtomIndex);

					bool hasFeasibleCoordinationComposition(const OriginalAtomIndex) const;
					bool hasInfeasibleChemicalBonds(const OriginalAtomIndex, const double errorRate) const;

					ChemicalComposition getCoordinationComposition(const OriginalAtomIndex) const;
					ChemicalComposition getCovalentCoordinationComposition(const OriginalAtomIndex) const;
					ChemicalComposition getIonicCoordinationComposition(const OriginalAtomIndex) const;
					ChemicalComposition toChemicalComposition(const std::vector<TranslatedAtomIndex>& translatedIndices) const;
					ChemicalComposition toCoordinationComposition(const std::vector<std::pair<double, OriginalAtomIndex>>&, const std::vector<std::pair<double, TranslatedAtomIndex>>&) const;

					std::vector<std::pair<double, OriginalAtomIndex>> getOrderedIonicBondedOriginalIndices(const OriginalAtomIndex centralIndex) const;
					std::vector<std::pair<double, TranslatedAtomIndex>> getOrderedIonicBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const;
					std::vector<std::pair<double, OriginalAtomIndex>> getOrderedCovalentBondedOriginalIndices(const OriginalAtomIndex centralIndex) const;
					std::vector<std::pair<double, TranslatedAtomIndex>> getOrderedCovalentBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const;

					std::vector<std::pair<double, OriginalAtomIndex>> getOrderedChemicalBondedOriginalIndices(const OriginalAtomIndex centralIndex) const;
					std::vector<std::pair<double, TranslatedAtomIndex>> getOrderedChemicalBondedTranslatedIndices(const OriginalAtomIndex centralIndex) const;
					std::vector<std::pair<double, OriginalAtomIndex>> getOrderedChemicalBondedOriginalIndices(const OriginalAtomIndex centralIndex, const AtomicNumber) const;
					std::vector<std::pair<double, TranslatedAtomIndex>> getOrderedChemicalBondedTranslatedIndices(const OriginalAtomIndex centralIndex, const AtomicNumber) const;


					std::vector<OriginalAtomIndex> getCovalentExcludedOriginalIndices(const OriginalAtomIndex) const;
					std::vector<TranslatedAtomIndex> getCovalentExcludedTranslatedIndices(const OriginalAtomIndex) const;
					std::vector<OriginalAtomIndex> getIonicExcludedOriginalIndices(const OriginalAtomIndex) const;
					std::vector<TranslatedAtomIndex> getIonicExcludedTranslatedIndices(const OriginalAtomIndex) const;

				// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

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
// Protected methods

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::clearChemicalBonds(const OriginalAtomIndex originalAtomIndex)
{
	const ConstrainingAtom& centralAtom = atoms()[originalAtomIndex];


	for (const auto& coordinatedIndex : centralAtom.getCovalentBondedOriginalAtomIndices())
		eraseCovalentBond(originalAtomIndex, coordinatedIndex);

	for (const auto& coordinatedIndex : centralAtom.getIonicBondedOriginalAtomIndices())
		eraseIonicBond(originalAtomIndex, coordinatedIndex);

	for (const auto& coordinatedIndex : centralAtom.getCovalentBondedTranslatedAtomIndices())
		eraseCovalentBond(originalAtomIndex, coordinatedIndex);

	for (const auto& coordinatedIndex : centralAtom.getIonicBondedTranslatedAtomIndices())
		eraseIonicBond(originalAtomIndex, coordinatedIndex);
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::getCoordinationComposition(const OriginalAtomIndex centralAtomIndex) const
{
	ChemicalComposition coordinationComposition;
	{
		const ConstrainingAtom& constrainingAtom = atoms()[centralAtomIndex];


		for (const auto originalIndex : constrainingAtom.getCovalentBondedOriginalAtomIndices())
			coordinationComposition.add(atoms()[originalIndex].ionicAtomicNumber().atomicNumber());

		for (const auto& translatedIndex : constrainingAtom.getCovalentBondedTranslatedAtomIndices())
			coordinationComposition.add(atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber());

		for (const auto originalIndex : constrainingAtom.getIonicBondedOriginalAtomIndices())
			coordinationComposition.add(atoms()[originalIndex].ionicAtomicNumber().atomicNumber());

		for (const auto& translatedIndex : constrainingAtom.getIonicBondedTranslatedAtomIndices())
			coordinationComposition.add(atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber());
	}

	return coordinationComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::getCovalentCoordinationComposition(const OriginalAtomIndex centralAtomIndex) const
{
	ChemicalComposition coordinationComposition;
	{
		const ConstrainingAtom& constrainingAtom = atoms()[centralAtomIndex];


		for (const auto originalIndex : constrainingAtom.getCovalentBondedOriginalAtomIndices())
			coordinationComposition.add(atoms()[originalIndex].ionicAtomicNumber().atomicNumber());

		for (const auto& translatedIndex : constrainingAtom.getCovalentBondedTranslatedAtomIndices())
			coordinationComposition.add(atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber());
	}

	return coordinationComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::getIonicCoordinationComposition(const OriginalAtomIndex centralAtomIndex) const
{
	ChemicalComposition coordinationComposition;
	{
		const ConstrainingAtom& constrainingAtom = atoms()[centralAtomIndex];


		for (const auto originalIndex : constrainingAtom.getIonicBondedOriginalAtomIndices())
			coordinationComposition.add(atoms()[originalIndex].ionicAtomicNumber().atomicNumber());

		for (const auto& translatedIndex : constrainingAtom.getIonicBondedTranslatedAtomIndices())
			coordinationComposition.add(atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber());
	}

	return coordinationComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::toChemicalComposition(const std::vector<TranslatedAtomIndex>& translatedIndices) const
{
	ChemicalComposition chemicalComposition;
	{
		for (const auto& translatedIndex : translatedIndices)
			chemicalComposition.add(atoms()[translatedIndex.originalIndex()].ionicAtomicNumber().atomicNumber());
	}

	return chemicalComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Components::Internal::CoordinationPolyhedraRetriever::toCoordinationComposition(const std::vector<std::pair<double, OriginalAtomIndex>>& originalIndices, const std::vector<std::pair<double, TranslatedAtomIndex>>& translatedIndices) const
{
	ChemicalComposition coordinationComposition;
	{
		for (const auto distanceAndIndex : originalIndices)
			coordinationComposition.add(atoms()[distanceAndIndex.second].ionicAtomicNumber().atomicNumber());

		for (const auto& distanceAndIndex : translatedIndices)
			coordinationComposition.add(atoms()[distanceAndIndex.second.originalIndex()].ionicAtomicNumber().atomicNumber());
	}

	return coordinationComposition;
}

// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_COORDINATIONPOLYHEDRARETRIEVER_H
