#include "FeasiblePolyhedraConnectionsDictionary.h"

#include <regex>

#include "CoordinationPolyhedraConnector.h"

using namespace MathematicalCrystalChemistry::CrystalModel::Constraints;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Methods

Internal::CoordinationPolyhedraConnectionParameters FeasiblePolyhedraConnectionsDictionary::s_coordinationPolyhedraConnectionParameters;
std::unordered_map<std::pair<FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber, FeasiblePolyhedraConnectionsDictionary::IonicAtomicNumber>, std::vector<FeasiblePolyhedraConnectionsDictionary::BasicChemicalComposition>, FeasiblePolyhedraConnectionsDictionary::KeyHasher, FeasiblePolyhedraConnectionsDictionary::KeyEqual> FeasiblePolyhedraConnectionsDictionary::s_feasiblePolyhedraLinkingDictionary;



void FeasiblePolyhedraConnectionsDictionary::initialize() noexcept
{
	s_coordinationPolyhedraConnectionParameters.initialize();

	s_feasiblePolyhedraLinkingDictionary.clear();
}

void FeasiblePolyhedraConnectionsDictionary::initialize(const System::IO::StreamReader& streamReader)
{
	s_coordinationPolyhedraConnectionParameters.initialize(streamReader);
}

void FeasiblePolyhedraConnectionsDictionary::initialize(const ChemicalComposition& crystalComposition)
{
	initialize(toConstrainingChemicalComposition(crystalComposition));
}

void FeasiblePolyhedraConnectionsDictionary::initialize(const ConstrainingChemicalComposition& constrainingCrystalComposition)
{
	Internal::CoordinationPolyhedraConnector coordinationPolyhedraConnector;
	coordinationPolyhedraConnector.setCoordinationPolyhedraConnectionParameters(s_coordinationPolyhedraConnectionParameters);


	s_feasiblePolyhedraLinkingDictionary.clear();
	{
		for (const auto& pairAndBridgings : coordinationPolyhedraConnector.getFeasibleLinkingDictionary(constrainingCrystalComposition))
		{
			std::vector<BasicChemicalComposition> bridgingCompositions;
			{
				for (const auto& bridgingComposition : pairAndBridgings.second)
					bridgingCompositions.push_back(bridgingComposition);
			}

			s_feasiblePolyhedraLinkingDictionary.emplace(pairAndBridgings.first, std::move(bridgingCompositions));
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
