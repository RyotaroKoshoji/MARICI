#ifndef MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOM_H
#define MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOM_H

#include <vector>
#include <unordered_set>

#include "NumericalVector.h"
#include "AtomIndex.h"
#include "IonicAtomicNumber.h"

#include "ConstrainingAtomicSpecies.h"


namespace MathematicalCrystalChemistry
{
	namespace CrystalModel
	{
		namespace Components
		{
			class OptimalAtom;
			class SphericalAtom;


			class ConstrainingAtom final
			{
				using size_type = unsigned short;
				using IonicAtomicNumber = ChemToolkit::Generic::IonicAtomicNumber;
				using NumericalVector = MathToolkit::LinearAlgebra::NumericalVector<double, 3>;

				using OriginalAtomIndex = ChemToolkit::Crystallography::OriginalAtomIndex;
				using TranslatedAtomIndex = ChemToolkit::Crystallography::TranslatedAtomIndex;

				using CoordinationConstraints = MathematicalCrystalChemistry::CrystalModel::Constraints::CoordinationConstraints;
				using CovalentRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::CovalentRadius;
				using IonicRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::IonicRadius;
				using IonicRepulsionRadius = MathematicalCrystalChemistry::CrystalModel::Constraints::IonicRepulsionRadius;


				struct OriginalHasher
				{
					std::size_t operator()(const OriginalAtomIndex) const noexcept;
				};

				struct TranslatedHasher
				{
					std::size_t operator()(const TranslatedAtomIndex&) const noexcept;
				};

// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Constructors, destructor, and operators

			public:
				ConstrainingAtom() noexcept;
				explicit ConstrainingAtom(const ConstrainingAtomicSpecies&) noexcept;
				ConstrainingAtom(const ConstrainingAtomicSpecies&, const NumericalVector&) noexcept;
				explicit ConstrainingAtom(const IonicAtomicNumber&);
				ConstrainingAtom(const IonicAtomicNumber&, const NumericalVector&);

				ConstrainingAtom(const IonicAtomicNumber&, const CoordinationConstraints&, const SphericalAtom&) noexcept;
				explicit ConstrainingAtom(const OptimalAtom&);

				~ConstrainingAtom() = default;

				ConstrainingAtom(const ConstrainingAtom&) = default;
				ConstrainingAtom(ConstrainingAtom&&) noexcept = default;
				ConstrainingAtom& operator=(const ConstrainingAtom&) = default;
				ConstrainingAtom& operator=(ConstrainingAtom&&) noexcept = default;

			// Constructors, destructor, and operators
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Property

				const IonicAtomicNumber& ionicAtomicNumber() const noexcept;
				const NumericalVector& cartesianCoordinate() const noexcept;
				NumericalVector& cartesianCoordinate() noexcept;

				const CoordinationConstraints& coordinationConstraints() const noexcept;
				const CovalentRadius& covalentRadius() const noexcept;
				const IonicRadius& ionicRadius() const noexcept;
				const IonicRepulsionRadius& ionicRepulsionRadius() const noexcept;

				std::vector<OriginalAtomIndex> getCovalentBondedOriginalAtomIndices() const noexcept;
				std::vector<OriginalAtomIndex> getIonicBondedOriginalAtomIndices() const noexcept;
				std::vector<OriginalAtomIndex> getIonicRepulsedOriginalAtomIndices() const noexcept;
				std::vector<TranslatedAtomIndex> getCovalentBondedTranslatedAtomIndices() const noexcept;
				std::vector<TranslatedAtomIndex> getIonicBondedTranslatedAtomIndices() const noexcept;
				std::vector<TranslatedAtomIndex> getIonicRepulsedTranslatedAtomIndices() const noexcept;

			// Property
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
			// Methods

				size_type getCoordinationNumber() const noexcept;
				size_type getCovalentCoordinationNumber() const noexcept;
				size_type getIonicCoordinationNumber() const noexcept;

				bool hasCovalentBondWith(const OriginalAtomIndex) const noexcept;
				bool hasCovalentBondWith(const TranslatedAtomIndex&) const noexcept;
				bool hasIonicBondWith(const OriginalAtomIndex) const noexcept;
				bool hasIonicBondWith(const TranslatedAtomIndex&) const noexcept;
				bool hasIonicRepulsionWith(const OriginalAtomIndex) const noexcept;
				bool hasIonicRepulsionWith(const TranslatedAtomIndex&) const noexcept;

				bool createCovalentBondWith(const OriginalAtomIndex) noexcept;
				bool createCovalentBondWith(const TranslatedAtomIndex&) noexcept;
				bool createIonicBondWith(const OriginalAtomIndex) noexcept;
				bool createIonicBondWith(const TranslatedAtomIndex&) noexcept;
				bool createIonicRepulsionWith(const OriginalAtomIndex) noexcept;
				bool createIonicRepulsionWith(const TranslatedAtomIndex&) noexcept;

				bool eraseCovalentBond(const OriginalAtomIndex) noexcept;
				bool eraseCovalentBond(const TranslatedAtomIndex&) noexcept;
				bool eraseIonicBond(const OriginalAtomIndex) noexcept;
				bool eraseIonicBond(const TranslatedAtomIndex&) noexcept;
				bool eraseIonicRepulsion(const OriginalAtomIndex) noexcept;
				bool eraseIonicRepulsion(const TranslatedAtomIndex&) noexcept;

				void clearCovalentBonds() noexcept;
				void clearIonicBonds() noexcept;
				void clearIonicRepulsions() noexcept;

			// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************

			private:
				IonicAtomicNumber _ionicAtomicNumber;
				NumericalVector _cartesianCoordinate;

				CoordinationConstraints _coordinationConstraints;
				CovalentRadius _covalentRadius;
				IonicRadius _ionicRadius;
				IonicRepulsionRadius _ionicRepulsionRadius;

				std::unordered_set<OriginalAtomIndex, OriginalHasher> _covalentBondedOriginalAtomIndices;
				std::unordered_set<OriginalAtomIndex, OriginalHasher> _ionicBondedOriginalAtomIndices;
				std::unordered_set<OriginalAtomIndex, OriginalHasher> _ionicRepulsedOriginalAtomIndices;
				std::unordered_set<TranslatedAtomIndex, TranslatedHasher> _covalentBondedTranslatedAtomIndices;
				std::unordered_set<TranslatedAtomIndex, TranslatedHasher> _ionicBondedTranslatedAtomIndices;
				std::unordered_set<TranslatedAtomIndex, TranslatedHasher> _ionicRepulsedTranslatedAtomIndices;
			};



			inline std::size_t ConstrainingAtom::OriginalHasher::operator()(const OriginalAtomIndex oai) const noexcept
			{
				return std::hash<OriginalAtomIndex>{}(oai);
			}

			inline std::size_t ConstrainingAtom::TranslatedHasher::operator()(const TranslatedAtomIndex& tai) const noexcept
			{
				return tai.getHashCode();
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

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::IonicAtomicNumber& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::ionicAtomicNumber() const noexcept
{
	return _ionicAtomicNumber;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::NumericalVector& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::cartesianCoordinate() const noexcept
{
	return _cartesianCoordinate;
}

inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::NumericalVector& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::cartesianCoordinate() noexcept
{
	return _cartesianCoordinate;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::CoordinationConstraints& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::coordinationConstraints() const noexcept
{
	return _coordinationConstraints;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::CovalentRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::covalentRadius() const noexcept
{
	return _covalentRadius;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::IonicRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::ionicRadius() const noexcept
{
	return _ionicRadius;
}

inline const MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::IonicRepulsionRadius& MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::ionicRepulsionRadius() const noexcept
{
	return _ionicRepulsionRadius;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::OriginalAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getCovalentBondedOriginalAtomIndices() const noexcept
{
	std::vector<OriginalAtomIndex> indices;
	{
		for (const auto& index : _covalentBondedOriginalAtomIndices)
			indices.push_back(index);
	}

	return indices;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::OriginalAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getIonicBondedOriginalAtomIndices() const noexcept
{
	std::vector<OriginalAtomIndex> indices;
	{
		for (const auto& index : _ionicBondedOriginalAtomIndices)
			indices.push_back(index);
	}

	return indices;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::OriginalAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getIonicRepulsedOriginalAtomIndices() const noexcept
{
	std::vector<OriginalAtomIndex> indices;
	{
		for (const auto& index : _ionicRepulsedOriginalAtomIndices)
			indices.push_back(index);
	}

	return indices;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::TranslatedAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getCovalentBondedTranslatedAtomIndices() const noexcept
{
	std::vector<TranslatedAtomIndex> indices;
	{
		for (const auto& index : _covalentBondedTranslatedAtomIndices)
			indices.push_back(index);
	}

	return indices;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::TranslatedAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getIonicBondedTranslatedAtomIndices() const noexcept
{
	std::vector<TranslatedAtomIndex> indices;
	{
		for (const auto& index : _ionicBondedTranslatedAtomIndices)
			indices.push_back(index);
	}

	return indices;
}

inline std::vector<MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::TranslatedAtomIndex> MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getIonicRepulsedTranslatedAtomIndices() const noexcept
{
	std::vector<TranslatedAtomIndex> indices;
	{
		for (const auto& index : _ionicRepulsedTranslatedAtomIndices)
			indices.push_back(index);
	}

	return indices;
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

inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getCoordinationNumber() const noexcept
{
	return (getCovalentCoordinationNumber() + getIonicCoordinationNumber());
}

inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getCovalentCoordinationNumber() const noexcept
{
	return static_cast<size_type>(_covalentBondedOriginalAtomIndices.size() + _covalentBondedTranslatedAtomIndices.size());
}

inline MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::size_type MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::getIonicCoordinationNumber() const noexcept
{
	return static_cast<size_type>(_ionicBondedOriginalAtomIndices.size() + _ionicBondedTranslatedAtomIndices.size());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasCovalentBondWith(const OriginalAtomIndex index) const noexcept
{
	return !(_covalentBondedOriginalAtomIndices.find(index) == _covalentBondedOriginalAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasCovalentBondWith(const TranslatedAtomIndex& index) const noexcept
{
	return !(_covalentBondedTranslatedAtomIndices.find(index) == _covalentBondedTranslatedAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasIonicBondWith(const OriginalAtomIndex index) const noexcept
{
	return !(_ionicBondedOriginalAtomIndices.find(index) == _ionicBondedOriginalAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasIonicBondWith(const TranslatedAtomIndex& index) const noexcept
{
	return !(_ionicBondedTranslatedAtomIndices.find(index) == _ionicBondedTranslatedAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasIonicRepulsionWith(const OriginalAtomIndex index) const noexcept
{
	return !(_ionicRepulsedOriginalAtomIndices.find(index) == _ionicRepulsedOriginalAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::hasIonicRepulsionWith(const TranslatedAtomIndex& index) const noexcept
{
	return !(_ionicRepulsedTranslatedAtomIndices.find(index) == _ionicRepulsedTranslatedAtomIndices.end());
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createCovalentBondWith(const OriginalAtomIndex index) noexcept
{
	return _covalentBondedOriginalAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createCovalentBondWith(const TranslatedAtomIndex& index) noexcept
{
	return _covalentBondedTranslatedAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createIonicBondWith(const OriginalAtomIndex index) noexcept
{
	return _ionicBondedOriginalAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createIonicBondWith(const TranslatedAtomIndex& index) noexcept
{
	return _ionicBondedTranslatedAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createIonicRepulsionWith(const OriginalAtomIndex index) noexcept
{
	return _ionicRepulsedOriginalAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::createIonicRepulsionWith(const TranslatedAtomIndex& index) noexcept
{
	return _ionicRepulsedTranslatedAtomIndices.emplace(index).second;
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseCovalentBond(const OriginalAtomIndex index) noexcept
{
	return static_cast<bool>(_covalentBondedOriginalAtomIndices.erase(index));
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseCovalentBond(const TranslatedAtomIndex& index) noexcept
{
	return static_cast<bool>(_covalentBondedTranslatedAtomIndices.erase(index));
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseIonicBond(const OriginalAtomIndex index) noexcept
{
	return static_cast<bool>(_ionicBondedOriginalAtomIndices.erase(index));
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseIonicBond(const TranslatedAtomIndex& index) noexcept
{
	return static_cast<bool>(_ionicBondedTranslatedAtomIndices.erase(index));
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseIonicRepulsion(const OriginalAtomIndex index) noexcept
{
	return static_cast<bool>(_ionicRepulsedOriginalAtomIndices.erase(index));
}

inline bool MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::eraseIonicRepulsion(const TranslatedAtomIndex& index) noexcept
{
	return static_cast<bool>(_ionicRepulsedTranslatedAtomIndices.erase(index));
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::clearCovalentBonds() noexcept
{
	_covalentBondedOriginalAtomIndices.clear();
	_covalentBondedTranslatedAtomIndices.clear();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::clearIonicBonds() noexcept
{
	_ionicBondedOriginalAtomIndices.clear();
	_ionicBondedTranslatedAtomIndices.clear();
}

inline void MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingAtom::clearIonicRepulsions() noexcept
{
	_ionicRepulsedOriginalAtomIndices.clear();
	_ionicRepulsedTranslatedAtomIndices.clear();
}

// Methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************


#endif // !MATHEMATICALCRYSTALCHEMISTRY_CRYSTALMODEL_COMPONENTS_CONSTRAININGATOM_H
