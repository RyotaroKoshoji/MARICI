#include "CoordinationConstraintsDictionary.h"

#include <regex>

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

std::unordered_map<CoordinationConstraintsDictionary::IonicAtomicNumber, std::vector<CoordinationConstraintsDictionary::ChemicalComposition>, CoordinationConstraintsDictionary::KeyHasher, CoordinationConstraintsDictionary::KeyEqual> CoordinationConstraintsDictionary::s_feasibleCompositionDictionary;
std::unordered_map<CoordinationConstraintsDictionary::IonicAtomicNumber, std::vector<CoordinationConstraintsDictionary::size_type>, CoordinationConstraintsDictionary::KeyHasher, CoordinationConstraintsDictionary::KeyEqual> CoordinationConstraintsDictionary::s_feasibleCovalentCoordinationNumbers;
std::unordered_map<CoordinationConstraintsDictionary::IonicAtomicNumber, std::vector<CoordinationConstraintsDictionary::size_type>, CoordinationConstraintsDictionary::KeyHasher, CoordinationConstraintsDictionary::KeyEqual> CoordinationConstraintsDictionary::s_feasibleIonicCoordinationNumbers;
std::unordered_map<CoordinationConstraintsDictionary::IonicAtomicNumber, std::vector<CoordinationConstraintsDictionary::ChemicalComposition>, CoordinationConstraintsDictionary::KeyHasher, CoordinationConstraintsDictionary::KeyEqual> CoordinationConstraintsDictionary::s_lowerBoundCompositionDictionary;



void CoordinationConstraintsDictionary::initialize() noexcept
{
	s_feasibleCompositionDictionary.clear();

	s_feasibleCovalentCoordinationNumbers.clear();
	s_feasibleIonicCoordinationNumbers.clear();
	s_lowerBoundCompositionDictionary.clear();
}

void CoordinationConstraintsDictionary::initialize(const System::IO::StreamReader& streamReader)
{
	System::IO::StreamReader coordinationCompositionReader = streamReader.getListBlock("&", "FEASIBLE_COORDINATION_COMPOSITIONS");
	{
		for (const auto& stringLine : coordinationCompositionReader.allLines())
		{
			std::vector<std::string> parameterTuple = System::IO::StreamReader::toParameterTuple(stringLine);


			IonicAtomicNumber ionicAtomicNumber = IonicAtomicNumber::toIonicAtomicNumber(parameterTuple.at(0));
			std::vector<ChemicalComposition> coordinationCompositions;
			{
				for (size_type index = 1; index < parameterTuple.size(); ++index)
					coordinationCompositions.push_back(toChemicalComposition(parameterTuple[index]));
			}

			s_feasibleCompositionDictionary.emplace(ionicAtomicNumber, std::move(coordinationCompositions));
		}
	}


	System::IO::StreamReader covalentCoordinationNumberReader = streamReader.getListBlock("&", "FEASIBLE_COVALENT_COORDINATION_NUMBERS");
	{
		for (const auto& stringLine : covalentCoordinationNumberReader.allLines())
		{
			std::vector<std::string> parameterTuple = System::IO::StreamReader::toParameterTuple(stringLine);


			IonicAtomicNumber ionicAtomicNumber = IonicAtomicNumber::toIonicAtomicNumber(parameterTuple.at(0));
			std::vector<size_type> coordinationNumbers;
			{
				for (size_type index = 1; index < parameterTuple.size(); ++index)
					coordinationNumbers.push_back(toCoordinationNumber(parameterTuple[index]));
			}

			s_feasibleCovalentCoordinationNumbers.emplace(ionicAtomicNumber, std::move(coordinationNumbers));
		}
	}

	
	System::IO::StreamReader ionicCoordinationNumberReader = streamReader.getListBlock("&", "FEASIBLE_IONIC_COORDINATION_NUMBERS");
	{
		for (const auto& stringLine : ionicCoordinationNumberReader.allLines())
		{
			std::vector<std::string> parameterTuple = System::IO::StreamReader::toParameterTuple(stringLine);


			IonicAtomicNumber ionicAtomicNumber = IonicAtomicNumber::toIonicAtomicNumber(parameterTuple.at(0));
			std::vector<size_type> coordinationNumbers;
			{
				for (size_type index = 1; index < parameterTuple.size(); ++index)
					coordinationNumbers.push_back(toCoordinationNumber(parameterTuple[index]));
			}

			s_feasibleIonicCoordinationNumbers.emplace(ionicAtomicNumber, std::move(coordinationNumbers));
		}
	}


	System::IO::StreamReader lowerBoundReader = streamReader.getListBlock("&", "LOWER_BOUND_COORDINATION_COMPOSITIONS");
	{
		for (const auto& stringLine : lowerBoundReader.allLines())
		{
			std::vector<std::string> parameterTuple = System::IO::StreamReader::toParameterTuple(stringLine);


			IonicAtomicNumber ionicAtomicNumber = IonicAtomicNumber::toIonicAtomicNumber(parameterTuple.at(0));
			std::vector<ChemicalComposition> lowerBoundCompositions;
			{
				for (size_type index = 1; index < parameterTuple.size(); ++index)
					lowerBoundCompositions.push_back(toChemicalComposition(parameterTuple[index]));
			}

			s_lowerBoundCompositionDictionary.emplace(ionicAtomicNumber, std::move(lowerBoundCompositions));
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
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Private methods

CoordinationConstraintsDictionary::size_type CoordinationConstraintsDictionary::toCoordinationNumber(const std::string& coordinationTexts)
{
	std::regex numPat{ "([\\d]+)" };
	std::smatch numSm;

	if (std::regex_match(coordinationTexts, numSm, numPat))
		return static_cast<size_type>(std::stoi(numSm.str(1)));
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ "ChemToolkit::Generic::ChemicalComposition::toCoordinationNumber", "Could not read coordination number." };
}

CoordinationConstraintsDictionary::ChemicalComposition CoordinationConstraintsDictionary::toChemicalComposition(const std::string& compositionTexts)
{
	ChemicalComposition chemicalComposition;
	{
		std::regex totalPat("([[:alnum:]]+[\\_]{0,1}[\\d]+[\\_]{0,1})*([[:alpha:]]+[\\_]{0,1}[\\d]+){1}");
		std::smatch totalSm;


		if (std::regex_match(compositionTexts, totalSm, totalPat))
		{
			std::regex patA("([[:alpha:]]+)([\\d]+)");
			std::regex patB("([[:alpha:]]+)[\\_]{1}([\\d]+)");


			for (std::sregex_iterator iter(compositionTexts.begin(), compositionTexts.end(), patA); iter != std::sregex_iterator{}; ++iter)
			{
				ChemToolkit::Generic::ElementSymbol elementSymbol{ iter->str(1) };
				size_type numAtoms = static_cast<size_type>(std::stoi(iter->str(2)));

				if (0 < numAtoms)
					chemicalComposition.add(elementSymbol.toAtomicNumber(), numAtoms);
				else
					throw System::ExceptionServices::ArgumentOutOfRangeException{ "ChemToolkit::Generic::ChemicalComposition::toChemicalComposition", "There is a zero composition." };
			}

			for (std::sregex_iterator iter(compositionTexts.begin(), compositionTexts.end(), patB); iter != std::sregex_iterator{}; ++iter)
			{
				ChemToolkit::Generic::ElementSymbol elementSymbol{ iter->str(1) };
				size_type numAtoms = static_cast<size_type>(std::stoi(iter->str(2)));


				if (0 < numAtoms)
					chemicalComposition.add(elementSymbol.toAtomicNumber(), numAtoms);
				else
					throw System::ExceptionServices::ArgumentOutOfRangeException{ "ChemToolkit::Generic::ChemicalComposition::toChemicalComposition", "There is a zero composition." };
			}
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "ChemToolkit::Generic::ChemicalComposition::toChemicalComposition", "Could not read chemical composition." };
	}

	return chemicalComposition;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
