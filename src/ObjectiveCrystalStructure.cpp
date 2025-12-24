#include "ObjectiveCrystalStructure.h"

#include "ArgumentOutOfRangeException.h"

#include "ConstrainingCrystalStructure.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

ObjectiveCrystalStructure::ObjectiveCrystalStructure() noexcept
	: CrystalStructure{}
	, _correspondingIonicAtomicNumbers{}
	, _correspondingCoordinationConstraints{}
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
	, _translatedCovalentBondedIndices{}
	, _translatedCovalentExcludedIndices{}
	, _translatedIonicBondedIndices{}
	, _translatedIonicExcludedIndices{}
	, _translatedIonicRepulsedIndices{}
{
}

ObjectiveCrystalStructure::ObjectiveCrystalStructure(const UnitCell& cell, const std::vector<SphericalAtom>& atoms, const std::vector<IonicAtomicNumber>& numbers, const std::vector<CoordinationConstraints>& coordinations) noexcept
	: CrystalStructure{ cell, atoms }
	, _correspondingIonicAtomicNumbers{ numbers }
	, _correspondingCoordinationConstraints{ coordinations }
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
	, _translatedCovalentBondedIndices{}
	, _translatedCovalentExcludedIndices{}
	, _translatedIonicBondedIndices{}
	, _translatedIonicExcludedIndices{}
	, _translatedIonicRepulsedIndices{}
{
}

ObjectiveCrystalStructure::ObjectiveCrystalStructure(const UnitCell& cell, std::vector<SphericalAtom>&& atoms, std::vector<IonicAtomicNumber>&& numbers, std::vector<CoordinationConstraints>&& coordinations) noexcept
	: CrystalStructure{ cell, std::move(atoms) }
	, _correspondingIonicAtomicNumbers{ std::move(numbers) }
	, _correspondingCoordinationConstraints{ std::move(coordinations) }
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
	, _translatedCovalentBondedIndices{}
	, _translatedCovalentExcludedIndices{}
	, _translatedIonicBondedIndices{}
	, _translatedIonicExcludedIndices{}
	, _translatedIonicRepulsedIndices{}
{
}

ObjectiveCrystalStructure::ObjectiveCrystalStructure(const ConstrainingCrystalStructure& structure)
	: CrystalStructure{ structure.unitCell() }
	, _correspondingIonicAtomicNumbers{}
	, _correspondingCoordinationConstraints{}
	, _covalentBondedIndices{}
	, _covalentExcludedIndices{}
	, _ionicBondedIndices{}
	, _ionicExcludedIndices{}
	, _ionicRepulsedIndices{}
	, _translatedCovalentBondedIndices{}
	, _translatedCovalentExcludedIndices{}
	, _translatedIonicBondedIndices{}
	, _translatedIonicExcludedIndices{}
	, _translatedIonicRepulsedIndices{}
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


	for (const auto& indices : structure.constrainingIndexPairs())
	{
		const ConstrainingAtom& constrainingAtom = structure.atoms()[indices.originalAtomIndex()];


		if (structure.isIonicAttractive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasIonicBondWith(indices.translatedAtomIndex()))
				_translatedIonicBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedIonicExcludedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
		}

		else if (structure.isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
				_translatedCovalentBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedIonicRepulsedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
		}

		else
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
				_translatedCovalentBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedCovalentExcludedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
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

bool ObjectiveCrystalStructure::isFeasible(const double feasibleErrorRate, const double exclusionRatio) const
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


	for (const auto& indices : _translatedCovalentBondedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().minimum() + atoms()[indices.translatedAtomIndex().originalIndex()].covalentRadius().minimum();
		double maximumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().maximum() + atoms()[indices.translatedAtomIndex().originalIndex()].covalentRadius().maximum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, maximumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _translatedCovalentExcludedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].covalentRadius().maximum() + atoms()[indices.translatedAtomIndex().originalIndex()].covalentRadius().maximum();
		minimumDistance *= exclusionRatio;

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _translatedIonicBondedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().minimum() + atoms()[indices.translatedAtomIndex().originalIndex()].ionicRadius().minimum();
		double maximumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().maximum() + atoms()[indices.translatedAtomIndex().originalIndex()].ionicRadius().maximum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, maximumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _translatedIonicExcludedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRadius().maximum() + atoms()[indices.translatedAtomIndex().originalIndex()].ionicRadius().maximum();
		minimumDistance *= exclusionRatio;

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}

	for (const auto& indices : _translatedIonicRepulsedIndices)
	{
		double minimumDistance = atoms()[indices.originalAtomIndex()].ionicRepulsionRadius().minimum() + atoms()[indices.translatedAtomIndex().originalIndex()].ionicRepulsionRadius().minimum();

		if (!(isFeasibleDistance(indices.originalAtomIndex(), indices.translatedAtomIndex(), minimumDistance, feasibleErrorRate)))
			return false;
	}

	return true;
}

void ObjectiveCrystalStructure::import(const ConstrainingCrystalStructure& structure)
{
	atoms().clear();
	_correspondingIonicAtomicNumbers.clear();
	_correspondingCoordinationConstraints.clear();

	_covalentBondedIndices.clear();
	_covalentExcludedIndices.clear();
	_ionicBondedIndices.clear();
	_ionicExcludedIndices.clear();
	_ionicRepulsedIndices.clear();
	_translatedCovalentBondedIndices.clear();
	_translatedCovalentExcludedIndices.clear();
	_translatedIonicBondedIndices.clear();
	_translatedIonicExcludedIndices.clear();
	_translatedIonicRepulsedIndices.clear();


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


	for (const auto& indices : structure.constrainingIndexPairs())
	{
		const ConstrainingAtom& constrainingAtom = structure.atoms()[indices.originalAtomIndex()];


		if (structure.isIonicAttractive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasIonicBondWith(indices.translatedAtomIndex()))
				_translatedIonicBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedIonicExcludedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
		}

		else if (structure.isIonicRepulsive(indices.originalAtomIndex(), indices.translatedAtomIndex().originalIndex()))
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
				_translatedCovalentBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedIonicRepulsedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
		}

		else
		{
			if (constrainingAtom.hasCovalentBondWith(indices.translatedAtomIndex()))
				_translatedCovalentBondedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
			else
				_translatedCovalentExcludedIndices.push_back(ConstrainerIndices<TranslatedAtomIndex>{ indices.originalAtomIndex(), indices.translatedAtomIndex() });
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
