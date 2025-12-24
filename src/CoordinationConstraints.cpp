#include "CoordinationConstraints.h"

#include <iterator>

#include "CoordinationConstraintsDictionary.h"

#include "InvalidFileException.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Constraints;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

CoordinationConstraints::CoordinationConstraints() noexcept
	: _feasibleCompositionDictionary{}
	, _lowerBoundCompositionDictionary{}
	, _feasibleCovalentCoordinationNumbers{}
	, _feasibleIonicCoordinationNumbers{}
	, _componentAtomicNumbers{}
	, _maxCoordinationNumber{ std::numeric_limits<size_type>::max() }
	, _maxConstrainedCovalentCoordinationNumber{ std::numeric_limits<size_type>::max() }
	, _maxConstrainedIonicCoordinationNumber{ std::numeric_limits<size_type>::max() }
{
}

CoordinationConstraints::CoordinationConstraints(const IonicAtomicNumber& ionicAtomicNumber)
	: _feasibleCompositionDictionary{}
	, _lowerBoundCompositionDictionary{}
	, _feasibleCovalentCoordinationNumbers{}
	, _feasibleIonicCoordinationNumbers{}
	, _componentAtomicNumbers{}
	, _maxCoordinationNumber{}
	, _maxConstrainedCovalentCoordinationNumber{ std::numeric_limits<size_type>::max() }
	, _maxConstrainedIonicCoordinationNumber{ std::numeric_limits<size_type>::max() }
{
	if (CoordinationConstraintsDictionary::hasFeasibleCoordnationCompositions(ionicAtomicNumber))
	{
		if ((CoordinationConstraintsDictionary::hasLowerBoundCoordnationCompositions(ionicAtomicNumber)) || (CoordinationConstraintsDictionary::hasFeasibleCovalentCoordinationNumbers(ionicAtomicNumber)) || (CoordinationConstraintsDictionary::hasFeasibleIonicCoordinationNumbers(ionicAtomicNumber)))
			throw System::IO::InvalidFileException{ typeid(*this), "constructor", "Lower bound coordination compositions and/or Feasible covalent and/or ionic coordination numbers must be empty if there are feasible coordination compositions." };

		else
		{
			for (const auto& chemicalComposition : CoordinationConstraintsDictionary::getFeasibleCoordnationCompositions(ionicAtomicNumber))
			{
				for (const auto& numAndCount : chemicalComposition)
					_componentAtomicNumbers.emplace(numAndCount.first);

				_feasibleCompositionDictionary.emplace(chemicalComposition);
			}

			_maxCoordinationNumber = std::prev(_feasibleCompositionDictionary.end())->count();
		}
	}

	else
	{
		for (const auto& chemicalComposition : CoordinationConstraintsDictionary::getLowerBoundCoordnationCompositions(ionicAtomicNumber))
		{
			for (const auto& numAndCount : chemicalComposition)
				_componentAtomicNumbers.emplace(numAndCount.first);

			_lowerBoundCompositionDictionary.emplace(chemicalComposition);
		}

		for (const auto& coordinationNumber : CoordinationConstraintsDictionary::getFeasibleCovalentCoordinationNumbers(ionicAtomicNumber))
			_feasibleCovalentCoordinationNumbers.emplace(coordinationNumber);

		for (const auto& coordinationNumber : CoordinationConstraintsDictionary::getFeasibleIonicCoordinationNumbers(ionicAtomicNumber))
			_feasibleIonicCoordinationNumbers.emplace(coordinationNumber);


		{
			if (_feasibleCovalentCoordinationNumbers.empty())
				_maxConstrainedCovalentCoordinationNumber = std::numeric_limits<size_type>::max();
			else
				_maxConstrainedCovalentCoordinationNumber = *(std::prev(_feasibleCovalentCoordinationNumbers.end()));
		}
		{
			if (_feasibleIonicCoordinationNumbers.empty())
				_maxConstrainedIonicCoordinationNumber = std::numeric_limits<size_type>::max();
			else
				_maxConstrainedIonicCoordinationNumber = *(std::prev(_feasibleIonicCoordinationNumbers.end()));
		}
		_maxCoordinationNumber = _maxConstrainedCovalentCoordinationNumber + _maxConstrainedIonicCoordinationNumber;


		if (!(_lowerBoundCompositionDictionary.empty()))
		{
			if (_maxCoordinationNumber < _lowerBoundCompositionDictionary.begin()->count())
				throw System::IO::InvalidFileException{ typeid(*this), "constructor", "Minimum number of atoms in lower bound coordination compositions is more than maximum coordination number." };
		}
	}
}

CoordinationConstraints::CoordinationConstraints(const AtomicNumber number, const charge_type formalCharge)
	: CoordinationConstraints{ IonicAtomicNumber{ number, formalCharge } }
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
// Methods

CoordinationConstraints::ChemicalComposition CoordinationConstraints::getClosestFeasibleComposition(const ChemicalComposition& argumentComposition) const
{
	if (_feasibleCompositionDictionary.empty())
		return argumentComposition;

	else
	{
		ChemicalComposition closestComposition;
		{
			for (const auto& feasibleComposition : _feasibleCompositionDictionary)
			{
				if (argumentComposition < feasibleComposition)
					break;

				else if (feasibleComposition == argumentComposition)
					return feasibleComposition;

				else
				{
					bool isLowerFeasibleComposition = true;
					{
						for (const auto& numAndCount : feasibleComposition)
						{
							if (argumentComposition.count(numAndCount.first) < numAndCount.second)
							{
								isLowerFeasibleComposition = false;
								break;
							}
						}
					}

					if (isLowerFeasibleComposition)
						closestComposition = feasibleComposition;
				}
			}
		}

		return closestComposition;
	}
}

CoordinationConstraints::ChemicalComposition CoordinationConstraints::getClosestLowerBoundComposition(const ChemicalComposition& covalentComposition, const ChemicalComposition& ionicComposition) const
{
	if (_lowerBoundCompositionDictionary.empty())
	{
		ChemicalComposition closestLowerBoundComposition = covalentComposition;
		closestLowerBoundComposition.join(ionicComposition);

		return closestLowerBoundComposition;
	}
		
	else
	{
		ChemicalComposition closestLowerBoundComposition;
		{
			size_type closestCovalentCoordinationNumber = getClosestFeasibleCovalentCoordinationNumber(covalentComposition.count());
			size_type closestIonicCoordinationNumber = getClosestFeasibleIonicCoordinationNumber(ionicComposition.count());


			for (const auto& lowerBoundComposition : _lowerBoundCompositionDictionary)
			{
				if ((closestCovalentCoordinationNumber + closestIonicCoordinationNumber) < lowerBoundComposition.count())
					break;

				else
				{
					size_type covalentBonding = 0;
					size_type ionicBonding = 0;
					bool isLowerFeasibleComposition = true;
					{
						for (const auto& numAndCount : lowerBoundComposition)
						{
							if (numAndCount.second <= covalentComposition.count(numAndCount.first))
								++covalentBonding;

							else if (numAndCount.second <= ionicComposition.count(numAndCount.first))
								++ionicBonding;

							else
							{
								isLowerFeasibleComposition = false;
								break;
							}
						}
					}

					if (isLowerFeasibleComposition && (covalentBonding <= closestCovalentCoordinationNumber) && (ionicBonding < closestIonicCoordinationNumber))
						closestLowerBoundComposition = lowerBoundComposition;
				}
			}
		}

		return closestLowerBoundComposition;
	}
}

CoordinationConstraints::size_type CoordinationConstraints::getClosestFeasibleCovalentCoordinationNumber(const size_type argumentCoordinationNumber) const
{
	if (_feasibleCovalentCoordinationNumbers.empty())
		return argumentCoordinationNumber;

	else
	{
		size_type closestCoordinationNumber = 0;

		for (const auto num : _feasibleCovalentCoordinationNumbers)
		{
			if (argumentCoordinationNumber < num)
				break;
			else if (argumentCoordinationNumber == num)
				return num;
			else
				closestCoordinationNumber = num;
		}

		return closestCoordinationNumber;
	}
}

CoordinationConstraints::size_type CoordinationConstraints::getClosestFeasibleIonicCoordinationNumber(const size_type argumentCoordinationNumber) const
{
	if (_feasibleIonicCoordinationNumbers.empty())
		return argumentCoordinationNumber;

	else
	{
		size_type closestCoordinationNumber = 0;

		for (const auto num : _feasibleIonicCoordinationNumbers)
		{
			if (argumentCoordinationNumber < num)
				break;
			else if (argumentCoordinationNumber == num)
				return num;
			else
				closestCoordinationNumber = num;
		}

		return closestCoordinationNumber;
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
