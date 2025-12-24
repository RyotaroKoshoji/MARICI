#include "AtomicRadiusDictionary.h"

#include <regex>

#include "LengthCasting.h"

#include "InvalidFileException.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Constraints;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Methods

std::unordered_map<AtomicRadiusDictionary::IonicAtomicNumber, double, AtomicRadiusDictionary::KeyHasher, AtomicRadiusDictionary::KeyEqual> AtomicRadiusDictionary::s_minCovalentRadiusDictionary;
std::unordered_map<AtomicRadiusDictionary::IonicAtomicNumber, double, AtomicRadiusDictionary::KeyHasher, AtomicRadiusDictionary::KeyEqual> AtomicRadiusDictionary::s_maxCovalentRadiusDictionary;
std::unordered_map<AtomicRadiusDictionary::IonicAtomicNumber, double, AtomicRadiusDictionary::KeyHasher, AtomicRadiusDictionary::KeyEqual> AtomicRadiusDictionary::s_minIonicRadiusDictionary;
std::unordered_map<AtomicRadiusDictionary::IonicAtomicNumber, double, AtomicRadiusDictionary::KeyHasher, AtomicRadiusDictionary::KeyEqual> AtomicRadiusDictionary::s_maxIonicRadiusDictionary;
std::unordered_map<AtomicRadiusDictionary::IonicAtomicNumber, double, AtomicRadiusDictionary::KeyHasher, AtomicRadiusDictionary::KeyEqual> AtomicRadiusDictionary::s_minIonicRepulsionRadiusDictionary;



void AtomicRadiusDictionary::initialize() noexcept
{
	s_minCovalentRadiusDictionary.clear();
	s_maxCovalentRadiusDictionary.clear();

	s_minIonicRadiusDictionary.clear();
	s_maxIonicRadiusDictionary.clear();

	s_minIonicRepulsionRadiusDictionary.clear();
}

void AtomicRadiusDictionary::initialize(const System::IO::StreamReader& streamReader)
{
	std::regex pat;
	std::smatch sm;
	{
		std::string patTexts{ "[\\s\\t]*([\\+\\-[:alnum:]]+)[\\s\\t]+" };
		{
			for (size_type index = 0; index < 4; ++index)
				patTexts += "([\\.\\d]+)[\\s\\t]+";

			patTexts += "([\\.\\d]+)[\\s\\t]*";
		}

		pat = std::regex{ patTexts };
	}


	for (const auto& stringLine : streamReader.allLines())
	{
		if (std::regex_match(stringLine, sm, pat))
		{
			IonicAtomicNumber ionicAtomicNumber = IonicAtomicNumber::toIonicAtomicNumber(sm.str(1));

			using namespace MathToolkit::UnitConversion::LengthCasting;
			double minCovalentRadius = cast<Unit::Angstrom, Unit::AtomicUnit>(System::IO::StreamReader::toFractionalValue(sm.str(2)));
			double maxCovalentRadius = cast<Unit::Angstrom, Unit::AtomicUnit>(System::IO::StreamReader::toFractionalValue(sm.str(3)));
			double minIonicRadius = cast<Unit::Angstrom, Unit::AtomicUnit>(System::IO::StreamReader::toFractionalValue(sm.str(4)));
			double maxIonicRadius = cast<Unit::Angstrom, Unit::AtomicUnit>(System::IO::StreamReader::toFractionalValue(sm.str(5)));
			double minIonicRepulsionRadius = cast<Unit::Angstrom, Unit::AtomicUnit>(System::IO::StreamReader::toFractionalValue(sm.str(6)));


			validateCovalentRadii(minCovalentRadius, maxCovalentRadius);
			validateIonicRadii(minIonicRadius, maxIonicRadius);
			validateIonicRepulsionRadii(minIonicRepulsionRadius);


			s_minCovalentRadiusDictionary.emplace(ionicAtomicNumber, minCovalentRadius);
			s_maxCovalentRadiusDictionary.emplace(ionicAtomicNumber, maxCovalentRadius);
			s_minIonicRadiusDictionary.emplace(ionicAtomicNumber, minIonicRadius);
			s_maxIonicRadiusDictionary.emplace(ionicAtomicNumber, maxIonicRadius);
			s_minIonicRepulsionRadiusDictionary.emplace(ionicAtomicNumber, minIonicRepulsionRadius);
		}

		else
			throw System::IO::InvalidFileException{  "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::initialize", "Could not read atomic radii." };
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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private methods

void AtomicRadiusDictionary::validateCovalentRadii(const double minimum, const double maximum)
{
	if(minimum < 0.0)
		throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateCovalentRadii", "Minimum covalent radius is less than zero." };

	else
	{
		if (maximum < 0.0)
			throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateCovalentRadii", "Maximum covalent radius is less than zero." };

		else
		{
			if (maximum < minimum)
				throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateCovalentRadii", "Maximum covalent radius is less than minimum covalent radius." };
		}
	}
}

void AtomicRadiusDictionary::validateIonicRadii(const double minimum, const double maximum)
{
	if (minimum < 0.0)
		throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateIonicRadii", "Minimum ionic radius is less than zero." };

	else
	{
		if (maximum < 0.0)
			throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateIonicRadii", "Maximum ionic radius is less than zero." };

		else
		{
			if (maximum < minimum)
				throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateIonicRadii", "Maximum ionic radius is less than minimum ionic radius." };
		}
	}
}

void AtomicRadiusDictionary::validateIonicRepulsionRadii(const double minimum)
{
	if (minimum < 0.0)
		throw System::IO::InvalidFileException{ "MathematicalCrystalChemistry::CrystalModel::Constraints::AtomicRadiusDictionary::validateIonicRepulsionRadii", "Minimum ionic repulsion radius is less than zero." };
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
