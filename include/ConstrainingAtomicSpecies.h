#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOMSPECIES_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOMSPECIES_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include "AtomicNumber.h"
#include "ElementSymbol.h"

#include "AtomicRadius.h"
#include "CoordinationConstraints.h"

#include "IonicAtomicNumber.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			class ConstrainingAtomicSpecies final
			{
				using CoordinationConstraints = MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints;
				using CovalentRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::CovalentRadius;
				using IonicRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::IonicRadius;
				using IonicRepulsionRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::IonicRepulsionRadius;


			public:
				using size_type = unsigned short;
				using charge_type = short;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;

				struct Hasher
				{
					std::size_t operator()(const ConstrainingAtomicSpecies&) const noexcept;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				ConstrainingAtomicSpecies() noexcept;
				explicit ConstrainingAtomicSpecies(const IonicAtomicNumber&);

				~ConstrainingAtomicSpecies() = default;

				ConstrainingAtomicSpecies(const ConstrainingAtomicSpecies&) = default;
				ConstrainingAtomicSpecies(ConstrainingAtomicSpecies&&) noexcept = default;
				ConstrainingAtomicSpecies& operator=(const ConstrainingAtomicSpecies&) = default;
				ConstrainingAtomicSpecies& operator=(ConstrainingAtomicSpecies&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const IonicAtomicNumber& ionicAtomicNumber() const noexcept;
				const CoordinationConstraints& coordinationConstraints() const noexcept;

				const CovalentRadius& covalentRadius() const noexcept;
				const IonicRadius& ionicRadius() const noexcept;
				const IonicRepulsionRadius& ionicRepulsionRadius() const noexcept;


				void setIonicAtomicNumber(const IonicAtomicNumber&) noexcept;
				void setCoordinationConstraints(const CoordinationConstraints&) noexcept;

				void setCovalentRadius(const CovalentRadius&) noexcept;
				void setIonicRadius(const IonicRadius&) noexcept;
				void setIonicRepulsionRadius(const IonicRepulsionRadius&) noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				std::size_t getHashCode() const;
				std::string toHashString() const;
				std::string toElementSymbol() const;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				IonicAtomicNumber _ionicAtomicNumber;
				CoordinationConstraints _coordinationConstraints;

				CovalentRadius _covalentRadius;
				IonicRadius _ionicRadius;
				IonicRepulsionRadius _ionicRepulsionRadius;
			};



			inline std::size_t ConstrainingAtomicSpecies::Hasher::operator()(const ConstrainingAtomicSpecies& cas) const noexcept
			{
				return cas.getHashCode();
			}

			inline bool operator<(const ConstrainingAtomicSpecies& casa, const ConstrainingAtomicSpecies& casb) noexcept
			{
				if (casa.ionicAtomicNumber() < casb.ionicAtomicNumber())
					return true;
				else
					return false;
			}

			inline bool operator>(const ConstrainingAtomicSpecies& casa, const ConstrainingAtomicSpecies& casb) noexcept
			{
				return (casb < casa);
			}

			inline bool operator==(const ConstrainingAtomicSpecies& casa, const ConstrainingAtomicSpecies& casb) noexcept
			{
				if ((casa.ionicAtomicNumber() == casb.ionicAtomicNumber()))
					return true;
				else
					return false;
			}

			inline bool operator!=(const ConstrainingAtomicSpecies& casa, const ConstrainingAtomicSpecies& casb) noexcept
			{
				return !(casa == casb);
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

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::IonicAtomicNumber& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::ionicAtomicNumber() const noexcept
{
	return _ionicAtomicNumber;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::CoordinationConstraints& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::coordinationConstraints() const noexcept
{
	return _coordinationConstraints;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::CovalentRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::covalentRadius() const noexcept
{
	return _covalentRadius;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::IonicRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::ionicRadius() const noexcept
{
	return _ionicRadius;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::IonicRepulsionRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::ionicRepulsionRadius() const noexcept
{
	return _ionicRepulsionRadius;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::setIonicAtomicNumber(const IonicAtomicNumber& number) noexcept
{
	_ionicAtomicNumber = number;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::setCoordinationConstraints(const CoordinationConstraints& constraints) noexcept
{
	_coordinationConstraints = constraints;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::setCovalentRadius(const CovalentRadius& radius) noexcept
{
	_covalentRadius = radius;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::setIonicRadius(const IonicRadius& radius) noexcept
{
	_ionicRadius = radius;
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::setIonicRepulsionRadius(const IonicRepulsionRadius& radius) noexcept
{
	_ionicRepulsionRadius = radius;
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

inline std::size_t MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::getHashCode() const
{
	return _ionicAtomicNumber.getHashCode();
}

inline std::string MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::toHashString() const
{
	return _ionicAtomicNumber.toHashString();
}

inline std::string MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtomicSpecies::toElementSymbol() const
{
	return _ionicAtomicNumber.toString();
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOMSPECIES_H
