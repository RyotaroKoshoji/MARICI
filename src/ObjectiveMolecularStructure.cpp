#include "ObjectiveMolecularStructure.h"

#include "ArgumentOutOfRangeException.h"

#include "ConstrainingMolecularStructure.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

ObjectiveMolecularStructure::ObjectiveMolecularStructure() noexcept
	: CrystalStructure{}
	, _correspondingIonicAtomicNumbers{}
	, _correspondingCoordinationConstraints{}
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
{
}

ObjectiveMolecularStructure::ObjectiveMolecularStructure(const UnitCell& cell, const std::vector<SphericalAtom>& atoms, const std::vector<IonicAtomicNumber>& numbers, const std::vector<CoordinationConstraints>& coordinations) noexcept
	: CrystalStructure{ cell, atoms }
	, _correspondingIonicAtomicNumbers{ numbers }
	, _correspondingCoordinationConstraints{ coordinations }
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
{
}

ObjectiveMolecularStructure::ObjectiveMolecularStructure(const UnitCell& cell, std::vector<SphericalAtom>&& atoms, std::vector<IonicAtomicNumber>&& numbers, std::vector<CoordinationConstraints>&& coordinations) noexcept
	: CrystalStructure{ cell, std::move(atoms) }
	, _correspondingIonicAtomicNumbers{ std::move(numbers) }
	, _correspondingCoordinationConstraints{ std::move(coordinations) }
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
{
}

ObjectiveMolecularStructure::ObjectiveMolecularStructure(const ConstrainingMolecularStructure& structure)
	: CrystalStructure{ structure.unitCell() }
	, _correspondingIonicAtomicNumbers{}
	, _correspondingCoordinationConstraints{}
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
{
	for (const auto& atom : structure.atoms())
	{
		atoms().push_back(SphericalAtom{ atom });
		_correspondingIonicAtomicNumbers.push_back(atom.ionicAtomicNumber());
		_correspondingCoordinationConstraints.push_back(atom.coordinationConstraints());
	}


	for (size_type originalIndex = 0; originalIndex < structure.atoms().size(); ++originalIndex)
	{
		const ConstrainingAtom& constrainingAtom = structure.atoms()[originalIndex];


		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < structure.atoms().size(); ++translatedOriginalIndex)
		{
			if (structure.isIonicAttractive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasIonicBondWith(translatedOriginalIndex))
					_ionicBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableIonicExclusionDistance(originalIndex, translatedOriginalIndex))
						_ionicExcludedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
				}
			}

			else if (structure.isIonicRepulsive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
					_covalentBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableIonicRepulsionDistance(originalIndex, translatedOriginalIndex))
						_ionicRepulsedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
				}
			}

			else
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
					_covalentBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableCovalentExclusionDistance(originalIndex, translatedOriginalIndex))
						_covalentExcludedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
				}
			}
		}
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

bool ObjectiveMolecularStructure::isFeasible(const double feasibleErrorRate, const double exclusionRatio) const
{
	for (const auto& indices : _covalentBondedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().minimum() + atoms()[indices.translatedAtomIndex()].covalentRadius().minimum();
		double maximumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().maximum() + atoms()[indices.translatedAtomIndex()].covalentRadius().maximum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, maximumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _covalentExcludedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().maximum() + atoms()[indices.translatedAtomIndex()].covalentRadius().maximum();
		minimumDistance *= exclusionRatio;

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _ionicBondedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().minimum() + atoms()[indices.translatedAtomIndex()].ionicRadius().minimum();
		double maximumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().maximum() + atoms()[indices.translatedAtomIndex()].ionicRadius().maximum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, maximumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _ionicExcludedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().maximum() + atoms()[indices.translatedAtomIndex()].ionicRadius().maximum();
		minimumDistance *= exclusionRatio;

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _ionicRepulsedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRepulsionRadius().minimum() + atoms()[indices.translatedAtomIndex()].ionicRepulsionRadius().minimum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}


	return true;
}

void ObjectiveMolecularStructure::import(const ConstrainingMolecularStructure& structure)
{
	atoms().clear();
	_correspondingIonicAtomicNumbers.clear();
	_correspondingCoordinationConstraints.clear();

	_covalentBondedIndices.clear();
	_covalentExcludedIndices.clear();
	_ionicBondedIndices.clear();
	_ionicExcludedIndices.clear();
	_ionicRepulsedIndices.clear();


	unitCell() = structure.unitCell();

	for (const auto& atom : structure.atoms())
	{
		atoms().push_back(SphericalAtom{ atom });
		_correspondingIonicAtomicNumbers.push_back(atom.ionicAtomicNumber());
		_correspondingCoordinationConstraints.push_back(atom.coordinationConstraints());
	}


	for (size_type originalIndex = 0; originalIndex < structure.atoms().size(); ++originalIndex)
	{
		const ConstrainingAtom& constrainingAtom = structure.atoms()[originalIndex];


		for (size_type translatedOriginalIndex = (1 + originalIndex); translatedOriginalIndex < structure.atoms().size(); ++translatedOriginalIndex)
		{
			if (structure.isIonicAttractive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasIonicBondWith(translatedOriginalIndex))
					_ionicBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableIonicExclusionDistance(originalIndex, translatedOriginalIndex))
						_ionicExcludedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
				}
			}

			else if (structure.isIonicRepulsive(originalIndex, translatedOriginalIndex))
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
					_covalentBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableIonicRepulsionDistance(originalIndex, translatedOriginalIndex))
						_ionicRepulsedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
				}
			}

			else
			{
				if (constrainingAtom.hasCovalentBondWith(translatedOriginalIndex))
					_covalentBondedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });

				else
				{
					if (structure.isConstrainableCovalentExclusionDistance(originalIndex, translatedOriginalIndex))
						_covalentExcludedIndices.push_back(ConstrainerIndices<OriginalAtomIndex>{ originalIndex, translatedOriginalIndex });
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
