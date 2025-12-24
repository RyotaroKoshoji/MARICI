#include "ConstrainingAtom.h"

#include "OptimalAtom.h"
#include "SphericalAtom.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Components;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

ConstrainingAtom::ConstrainingAtom() noexcept
	: _ionicAtomicNumber{}
	, _cartesianCoordinate{ 0.0,0.0,0.0 }
	, _coordinationConstraints{}
	, _covalentRadius{}
	, _ionicRadius{}
	, _ionicRepulsionRadius{}
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const ConstrainingAtomicSpecies& species) noexcept
	: _ionicAtomicNumber{ species.ionicAtomicNumber() }
	, _cartesianCoordinate{ 0.0,0.0,0.0 }
	, _coordinationConstraints{ species.coordinationConstraints() }
	, _covalentRadius{ species.covalentRadius() }
	, _ionicRadius{ species.ionicRadius() }
	, _ionicRepulsionRadius{ species.ionicRepulsionRadius() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const ConstrainingAtomicSpecies& species, const NumericalVector& coordinate) noexcept
	: _ionicAtomicNumber{ species.ionicAtomicNumber() }
	, _cartesianCoordinate{ coordinate }
	, _coordinationConstraints{ species.coordinationConstraints() }
	, _covalentRadius{ species.covalentRadius() }
	, _ionicRadius{ species.ionicRadius() }
	, _ionicRepulsionRadius{ species.ionicRepulsionRadius() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const IonicAtomicNumber& ian)
	: _ionicAtomicNumber{ ian }
	, _cartesianCoordinate{ 0.0,0.0,0.0 }
	, _coordinationConstraints{ ian.atomicNumber(), ian.formalCharge() }
	, _covalentRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _ionicRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _ionicRepulsionRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const IonicAtomicNumber& ian, const NumericalVector& coordinate)
	: _ionicAtomicNumber{ ian }
	, _cartesianCoordinate{ coordinate }
	, _coordinationConstraints{ ian.atomicNumber(), ian.formalCharge() }
	, _covalentRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _ionicRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _ionicRepulsionRadius{ ian.atomicNumber(), ian.formalCharge() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const IonicAtomicNumber& ian, const CoordinationConstraints& coordinations, const SphericalAtom& sphericalAtom) noexcept
	: _ionicAtomicNumber{ ian }
	, _cartesianCoordinate{ sphericalAtom.cartesianCoordinate() }
	, _coordinationConstraints{ coordinations }
	, _covalentRadius{ sphericalAtom.covalentRadius() }
	, _ionicRadius{ sphericalAtom.ionicRadius() }
	, _ionicRepulsionRadius{ sphericalAtom.ionicRepulsionRadius() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
{
}

ConstrainingAtom::ConstrainingAtom(const OptimalAtom& optimalAtom)
	: _ionicAtomicNumber{ optimalAtom.ionicAtomicNumber() }
	, _cartesianCoordinate{ optimalAtom.cartesianCoordinate() }
	, _coordinationConstraints{ optimalAtom.ionicAtomicNumber().atomicNumber(), optimalAtom.ionicAtomicNumber().formalCharge() }
	, _covalentRadius{ optimalAtom.ionicAtomicNumber().atomicNumber(), optimalAtom.ionicAtomicNumber().formalCharge() }
	, _ionicRadius{ optimalAtom.ionicAtomicNumber().atomicNumber(), optimalAtom.ionicAtomicNumber().formalCharge() }
	, _ionicRepulsionRadius{ optimalAtom.ionicAtomicNumber().atomicNumber(), optimalAtom.ionicAtomicNumber().formalCharge() }
	, _covalentBondedOriginalAtomIndices{}
	, _ionicBondedOriginalAtomIndices{}
	, _ionicRepulsedOriginalAtomIndices{}
	, _covalentBondedTranslatedAtomIndices{}
	, _ionicBondedTranslatedAtomIndices{}
	, _ionicRepulsedTranslatedAtomIndices{}
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
