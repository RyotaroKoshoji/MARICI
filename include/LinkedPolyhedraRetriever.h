#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_LINKEDPOLYHEDRARETRIEVER_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_LINKEDPOLYHEDRARETRIEVER_H

#include "IonicAtomicNumber.h"
#include "FeasiblePolyhedraConnections.h"

#include "CoordinationPolyhedraRetriever.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			namespace Internal
			{
				class LinkedPolyhedraRetriever :public CoordinationPolyhedraRetriever
				{
					using LinkingList = std::vector<std::pair<ConstrainerIndices<TranslatedAtomIndex>, std::vector<TranslatedAtomIndex>>>;
					using LinkingDictionary = std::unordered_map<ConstrainerIndices<TranslatedAtomIndex>, std::vector<TranslatedAtomIndex>, ConstrainerIndices<TranslatedAtomIndex>::Hasher>;

				protected:
					using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
					using FeasiblePolyhedraConnections = MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnections;

// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Constructors, destructor, and operators

				protected:
					LinkedPolyhedraRetriever() noexcept;
					explicit LinkedPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&) noexcept;
					LinkedPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&, const std::vector<ConstrainingAtom>&) noexcept;
					LinkedPolyhedraRetriever(const ChemToolkit::Crystallography::UnitCell&, std::vector<ConstrainingAtom>&&) noexcept;

				public:
					virtual ~LinkedPolyhedraRetriever() = default;

					LinkedPolyhedraRetriever(const LinkedPolyhedraRetriever&) = default;
					LinkedPolyhedraRetriever(LinkedPolyhedraRetriever&&) noexcept = default;
					LinkedPolyhedraRetriever& operator=(const LinkedPolyhedraRetriever&) = default;
					LinkedPolyhedraRetriever& operator=(LinkedPolyhedraRetriever&&) noexcept = default;

				// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Property

					const FeasiblePolyhedraConnections& feasiblePolyhedraConnections() const noexcept;
					FeasiblePolyhedraConnections& feasiblePolyhedraConnections() noexcept;

				// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Public methods

					LinkingList getCoordinationPolyhedraLinkings() const;

				// Public methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Protected methods

				protected:
					std::pair<IonicAtomicNumber, IonicAtomicNumber> toIonicAtomicNumberPair(const ConstrainerIndices<TranslatedAtomIndex>&) const;

				// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Private methods

				private:
					void addPolyhedraLinkings(const OriginalAtomIndex, const std::vector<OriginalAtomIndex>&, LinkingDictionary&) const;
					void addPolyhedraLinkings(const OriginalAtomIndex, const std::vector<TranslatedAtomIndex>&, LinkingDictionary&) const;
					void addPolyhedraLinkings(const OriginalAtomIndex, const std::vector<OriginalAtomIndex>&, const std::vector<TranslatedAtomIndex>&, LinkingDictionary&) const;

				// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Private utility

					void addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& formerBondedIndex, const TranslatedAtomIndex& latterBondedIndex, LinkingDictionary&) const;
					void addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex formerBondedIndex, const TranslatedAtomIndex& latterBondedIndex, LinkingDictionary&) const;
					void addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& latterBondedIndex, const OriginalAtomIndex formerBondedIndex, LinkingDictionary&) const;

				// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************

				private:
					FeasiblePolyhedraConnections _feasiblePolyhedraConnections;
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
// Property

inline const MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::FeasiblePolyhedraConnections& MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::feasiblePolyhedraConnections() const noexcept
{
	return _feasiblePolyhedraConnections;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::FeasiblePolyhedraConnections& MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::feasiblePolyhedraConnections() noexcept
{
	return _feasiblePolyhedraConnections;
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
// Private methods

inline std::pair<MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::IonicAtomicNumber, MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::IonicAtomicNumber> MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::toIonicAtomicNumberPair(const ConstrainerIndices<TranslatedAtomIndex>& constrainerIndices) const
{
	return std::make_pair(atoms()[constrainerIndices.originalAtomIndex()].ionicAtomicNumber(), atoms()[constrainerIndices.translatedAtomIndex().originalIndex()].ionicAtomicNumber());
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

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinkings(const OriginalAtomIndex originalAtomIndex, const std::vector<OriginalAtomIndex>& chemicalBondedOriginalAtomIndices, LinkingDictionary& linkingDictionary) const
{
	for (const auto& formerBondedIndex : chemicalBondedOriginalAtomIndices)
	{
		for (const auto& latterBondedIndex : chemicalBondedOriginalAtomIndices)
		{
			if (formerBondedIndex < latterBondedIndex)
			{
				ConstrainerIndices<TranslatedAtomIndex> constrainerIndices{ formerBondedIndex, TranslatedAtomIndex{ latterBondedIndex } };
				auto iter = linkingDictionary.find(constrainerIndices);

				if (iter == linkingDictionary.end())
					linkingDictionary.emplace(constrainerIndices, std::vector<TranslatedAtomIndex>{ TranslatedAtomIndex{ originalAtomIndex } });
				else
					iter->second.push_back(TranslatedAtomIndex{ originalAtomIndex });
			}
		}
	}
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinkings(const OriginalAtomIndex originalAtomIndex, const std::vector<TranslatedAtomIndex>& chemicalBondedTranslatedAtomIndices, LinkingDictionary& linkingDictionary) const
{
	for (const auto& formerBondedIndex : chemicalBondedTranslatedAtomIndices)
	{
		for (const auto& latterBondedIndex : chemicalBondedTranslatedAtomIndices)
		{
			if (formerBondedIndex.originalIndex() < latterBondedIndex.originalIndex())
				addPolyhedraLinking(originalAtomIndex, formerBondedIndex, latterBondedIndex, linkingDictionary);

			else
			{
				if ((formerBondedIndex.originalIndex() == latterBondedIndex.originalIndex()) && (formerBondedIndex.latticePoint() < latterBondedIndex.latticePoint()))
					addPolyhedraLinking(originalAtomIndex, formerBondedIndex, latterBondedIndex, linkingDictionary);
			}
		}
	}
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinkings(const OriginalAtomIndex originalAtomIndex, const std::vector<OriginalAtomIndex>& chemicalBondedOriginalAtomIndices, const std::vector<TranslatedAtomIndex>& chemicalBondedTranslatedAtomIndices, LinkingDictionary& linkingDictionary) const
{
	constexpr LatticePoint originalLatticePoint{ 0,0,0 };


	for (const auto& formerBondedIndex : chemicalBondedOriginalAtomIndices)
	{
		for (const auto& latterBondedIndex : chemicalBondedTranslatedAtomIndices)
		{
			if (formerBondedIndex < latterBondedIndex.originalIndex())
				addPolyhedraLinking(originalAtomIndex, formerBondedIndex, latterBondedIndex, linkingDictionary);

			else
			{
				if (latterBondedIndex.originalIndex() < formerBondedIndex)
					addPolyhedraLinking(originalAtomIndex, latterBondedIndex, formerBondedIndex, linkingDictionary);

				else
				{
					if (originalLatticePoint < latterBondedIndex.latticePoint())
						addPolyhedraLinking(originalAtomIndex, formerBondedIndex, latterBondedIndex, linkingDictionary);

					else if (latterBondedIndex.latticePoint() < originalLatticePoint)
						addPolyhedraLinking(originalAtomIndex, latterBondedIndex, formerBondedIndex, linkingDictionary);

					else
						throw System::ExceptionServices::InvalidOperationException{ typeid(*this), "addPolyhedraLinkings", "A translated bonded index is in the original cell." };
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

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& formerBondedIndex, const TranslatedAtomIndex& latterBondedIndex, LinkingDictionary& linkingDictionary) const
{
	TranslatedAtomIndex translatedLatterBondedIndex{ latterBondedIndex.originalIndex(), formerBondedIndex.getRelativeLatticePointTo(latterBondedIndex) };
	TranslatedAtomIndex translatedOriginalIndex{ originalAtomIndex, formerBondedIndex.latticePoint() };
	translatedOriginalIndex.reverseLatticePoint();


	ConstrainerIndices<TranslatedAtomIndex> constrainerIndices{ formerBondedIndex.originalIndex(), translatedLatterBondedIndex };
	auto iter = linkingDictionary.find(constrainerIndices);


	if (iter == linkingDictionary.end())
		linkingDictionary.emplace(constrainerIndices, std::vector<TranslatedAtomIndex>{ translatedOriginalIndex });
	else
		iter->second.push_back(std::move(translatedOriginalIndex));
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const OriginalAtomIndex formerBondedIndex, const TranslatedAtomIndex& latterBondedIndex, LinkingDictionary& linkingDictionary) const
{
	ConstrainerIndices<TranslatedAtomIndex> constrainerIndices{ formerBondedIndex, latterBondedIndex };
	auto iter = linkingDictionary.find(constrainerIndices);


	if (iter == linkingDictionary.end())
		linkingDictionary.emplace(constrainerIndices, std::vector<TranslatedAtomIndex>{ TranslatedAtomIndex{ originalAtomIndex } });
	else
		iter->second.push_back(TranslatedAtomIndex{ originalAtomIndex });
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::Internal::LinkedPolyhedraRetriever::addPolyhedraLinking(const OriginalAtomIndex originalAtomIndex, const TranslatedAtomIndex& latterBondedIndex, const OriginalAtomIndex formerBondedIndex, LinkingDictionary& linkingDictionary) const
{
	TranslatedAtomIndex translatedFormerBondedIndex{ formerBondedIndex, latterBondedIndex.latticePoint() };
	translatedFormerBondedIndex.reverseLatticePoint();

	TranslatedAtomIndex translatedOriginalIndex{ originalAtomIndex, latterBondedIndex.latticePoint() };
	translatedOriginalIndex.reverseLatticePoint();


	ConstrainerIndices<TranslatedAtomIndex> constrainerIndices{ latterBondedIndex.originalIndex(), translatedFormerBondedIndex };
	auto iter = linkingDictionary.find(constrainerIndices);


	if (iter == linkingDictionary.end())
		linkingDictionary.emplace(constrainerIndices, std::vector<TranslatedAtomIndex>{ translatedOriginalIndex });
	else
		iter->second.push_back(std::move(translatedOriginalIndex));
}

// Private utility
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_INTERNAL_LINKEDPOLYHEDRARETRIEVER_H
