#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_POLYHEDRALFEASIBILITYCHECKERBASE_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_POLYHEDRALFEASIBILITYCHECKERBASE_H

#include <unordered_map>

#include "AtomicNumber.h"
#include "ChemicalComposition.h"

#include "IonicAtomicNumber.h"

#include "ConstrainingAtomicSpecies.h"
#include "ConstrainingAtom.h"
#include "ConstrainingMolecularStructure.h"
#include "ObjectiveMolecularStructure.h"

#include "CoordinationPolyhedraConnectionParameters.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Constraints
		{
			namespace Internal
			{
				class PolyhedralFeasibilityCheckerBase
				{
				protected:
					using size_type = unsigned short;

					using AtomicNumber = ChemToolkit::Generic::AtomicNumber;
					using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
					using BasicChemicalComposition = ChemToolkit::Generic::ChemicalComposition<AtomicNumber>;
					using ChemicalComposition = ChemToolkit::Generic::ChemicalComposition<IonicAtomicNumber>;

					using ConstrainingAtomicSpecies = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies;
					using ConstrainingChemicalComposition = ChemToolkit::Generic::ChemicalComposition<ConstrainingAtomicSpecies>;

					using ConstrainingAtom = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom;
					using ConstrainingMolecularStructure = MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingMolecularStructure;
					using ObjectiveMolecularStructure = MathematicalCrystalChemistry::CrystalModel::Components::ObjectiveMolecularStructure;


// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Constructors, destructor, and operators

				protected:
					PolyhedralFeasibilityCheckerBase() noexcept = default;

				public:
					virtual ~PolyhedralFeasibilityCheckerBase() = default;

					PolyhedralFeasibilityCheckerBase(const PolyhedralFeasibilityCheckerBase&) = default;
					PolyhedralFeasibilityCheckerBase(PolyhedralFeasibilityCheckerBase&&) noexcept = default;
					PolyhedralFeasibilityCheckerBase& operator=(const PolyhedralFeasibilityCheckerBase&) = default;
					PolyhedralFeasibilityCheckerBase& operator=(PolyhedralFeasibilityCheckerBase&&) noexcept = default;

				// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
				// Protected methods

				protected:
					BasicChemicalComposition toBasicChemicalComposition(const ChemicalComposition&) const;
					BasicChemicalComposition toBasicChemicalComposition(const ConstrainingChemicalComposition&) const;

					ChemicalComposition toChemicalComposition(const BasicChemicalComposition&, const ChemicalComposition& crystalComposition) const;
					ChemicalComposition toChemicalComposition(const ConstrainingChemicalComposition&) const;

					ConstrainingChemicalComposition toConstrainingChemicalComposition(const ChemicalComposition&) const;
					ConstrainingChemicalComposition toConstrainingChemicalComposition(const BasicChemicalComposition&, const ConstrainingChemicalComposition& crystalComposition) const;

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

inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::BasicChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toBasicChemicalComposition(const ChemicalComposition& chemicalComposition) const
{
	BasicChemicalComposition basicChemicalComposition;
	{
		for (const auto& numAndCount : chemicalComposition)
			basicChemicalComposition.add(numAndCount.first.atomicNumber(), numAndCount.second);
	}

	return basicChemicalComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::BasicChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toBasicChemicalComposition(const ConstrainingChemicalComposition& constrainingChemicalComposition) const
{
	BasicChemicalComposition basicChemicalComposition;
	{
		for (const auto& speciesAndCount : constrainingChemicalComposition)
			basicChemicalComposition.add(speciesAndCount.first.ionicAtomicNumber().atomicNumber(), speciesAndCount.second);
	}

	return basicChemicalComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toChemicalComposition(const BasicChemicalComposition& basicChemicalComposition, const ChemicalComposition& crystalComposition) const
{
	ChemicalComposition chemicalComposition;
	{
		for (const auto& numAndCount : basicChemicalComposition)
		{
			for (const auto& crystalNumAndCount : crystalComposition)
			{
				if (crystalNumAndCount.first.atomicNumber() == numAndCount.first)
					chemicalComposition.add(crystalNumAndCount.first, numAndCount.second);
			}
		}
	}

	return chemicalComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::ChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toChemicalComposition(const ConstrainingChemicalComposition& constrainingChemicalComposition) const
{
	ChemicalComposition chemicalComposition;
	{
		for (const auto& numAndCount : constrainingChemicalComposition)
			chemicalComposition.add(numAndCount.first.ionicAtomicNumber(), numAndCount.second);
	}

	return chemicalComposition;
}

inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::ConstrainingChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toConstrainingChemicalComposition(const ChemicalComposition& chemicalComposition) const
{
	ConstrainingChemicalComposition constrainingChemicalComposition;
	{
		for (const auto& numAndCount : chemicalComposition)
			constrainingChemicalComposition.add(ConstrainingAtomicSpecies{ numAndCount.first }, numAndCount.second);
	}

	return constrainingChemicalComposition;
}


inline MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::ConstrainingChemicalComposition MathematicalCrystalChemistry::CrystalModel::Constraints::Internal::PolyhedralFeasibilityCheckerBase::toConstrainingChemicalComposition(const BasicChemicalComposition& basicChemicalComposition, const ConstrainingChemicalComposition& crystalComposition) const
{
	ConstrainingChemicalComposition constrainingChemicalComposition;
	{
		for (const auto& numAndCount : basicChemicalComposition)
		{
			for (const auto& crystalNumAndCount : crystalComposition)
			{
				if (crystalNumAndCount.first.ionicAtomicNumber().atomicNumber() == numAndCount.first)
					constrainingChemicalComposition.add(crystalNumAndCount.first, numAndCount.second);
			}
		}
	}

	return constrainingChemicalComposition;
}


// Protected methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_CONSTRAINTS_INTERNAL_POLYHEDRALFEASIBILITYCHECKERBASE_H
