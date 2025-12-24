#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALCRYSTALSTRUCTURE_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALCRYSTALSTRUCTURE_H

#include <algorithm>
#include <utility>
#include <vector>

#include "StreamReader.h"

#include "CrystallographicAtom.h"
#include "CrystallographicStructure.h"

#include "IonicAtomicNumber.h"

#include "OptimalAtom.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			class ConstrainingCrystalStructure;


			class OptimalCrystalStructure :public ChemToolkit::Crystallography::CrystallographicStructure<OptimalAtom>
			{
				using UnitCell = ChemToolkit::Crystallography::UnitCell;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				OptimalCrystalStructure() noexcept;
				OptimalCrystalStructure(const UnitCell&, const std::vector<OptimalAtom>&) noexcept;
				OptimalCrystalStructure(const UnitCell&, std::vector<OptimalAtom>&&) noexcept;

				explicit OptimalCrystalStructure(const CrystallographicStructure<ChemToolkit::Crystallography::CrystallographicAtom>&);
				explicit OptimalCrystalStructure(const ConstrainingCrystalStructure&);

				virtual ~OptimalCrystalStructure() = default;

				OptimalCrystalStructure(const OptimalCrystalStructure&) = default;
				OptimalCrystalStructure(OptimalCrystalStructure&&) noexcept = default;
				OptimalCrystalStructure& operator=(const OptimalCrystalStructure&) = default;
				OptimalCrystalStructure& operator=(OptimalCrystalStructure&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				void initialize() noexcept override;
				void sortAtoms();

				void import(const ConstrainingCrystalStructure&);

				std::string toAbridgedStructuralFingerprint() const;
				std::string toStructuralFingerprint() const;
				std::string toDetailedStructuralFingerprint() const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Private methods

			private:
				std::vector<std::pair<OptimalAtom, size_type>> getUniqueOptimalAtoms() const;

			// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			};


			
			inline bool operator<(const OptimalCrystalStructure& ocsa, const OptimalCrystalStructure& ocsb)
			{
				if (std::lexicographical_compare(ocsa.atoms().begin(), ocsa.atoms().end(), ocsb.atoms().begin(), ocsb.atoms().end()))
				{
					if (ocsa.spaceGroupNumber() < ocsb.spaceGroupNumber())
						return true;
					else
						return false;
				}

				else
					return false;
			}

			inline bool operator>(const OptimalCrystalStructure& ocsa, const OptimalCrystalStructure& ocsb)
			{
				return (ocsb < ocsa);
			}

			inline bool operator==(const OptimalCrystalStructure& ocsa, const OptimalCrystalStructure& ocsb)
			{
				if (ocsa.spaceGroupNumber() == ocsb.spaceGroupNumber())
				{
					auto iterA = ocsa.atoms().begin();
					auto iterB = ocsb.atoms().begin();
					auto lastA = ocsa.atoms().end();
					auto lastB = ocsb.atoms().end();
					{
						for (; (iterA != lastA && iterB != lastB); ++iterA, ++iterB)
						{
							if (*iterA != *iterB)
								return false;
						}
					}

					return true;
				}

				else
					return false;
			}

			inline bool operator!=(const OptimalCrystalStructure& ocsa, const OptimalCrystalStructure& ocsb)
			{
				return !(ocsb == ocsa);
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
// Methods

inline void MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure::sortAtoms()
{
	std::vector<OptimalAtom> atomicArrangement{ atoms() };
	std::sort(atomicArrangement.begin(), atomicArrangement.end());

	setAtoms(std::move(atomicArrangement));
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_OPTIMALCRYSTALSTRUCTURE_H
