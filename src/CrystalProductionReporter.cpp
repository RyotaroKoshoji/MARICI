#include "CrystalProductionReporter.h"

#include "Directory.h"
#include "DirectoryNotFoundException.h"
#include "FileStream.h"
#include "InvalidFileException.h"

#include "CrystallographicInformation.h"
#include "CifStreamWriter.h"

using namespace MathematicalCrystalChemistry::Design::Diagnostics;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

std::mutex CrystalProductionReporter::s_reportMutex;
std::string CrystalProductionReporter::s_fingerprintFilename{ "fingerprint.txt" };



CrystalProductionReporter::CrystalProductionReporter() noexcept
	: _crystalProductionReportParameters{}
	, _crystalDesignRecorder{}
{
}

CrystalProductionReporter::CrystalProductionReporter(const CrystalProductionReportParameters& parameters)
	: _crystalProductionReportParameters{ parameters }
	, _crystalDesignRecorder{ parameters.crystalDesignRecordParameters() }
{
	if (!(System::IO::Directory::exist(_crystalProductionReportParameters.designedCrystalsDirectoryPath())))
	{
		std::string message;
		{
			message += "Could not find the directory of \"";
			message += _crystalProductionReportParameters.designedCrystalsDirectoryPath().generic_string();
			message += "\".";
		}

		throw System::IO::DirectoryNotFoundException{ typeid(*this), "constructor", message };
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
// Property

void CrystalProductionReporter::setCrystalProductionReportParameters(const CrystalProductionReportParameters& parameters)
{
	_crystalProductionReportParameters = parameters;
	{
		if (!(System::IO::Directory::exist(_crystalProductionReportParameters.designedCrystalsDirectoryPath())))
		{
			std::string message;
			{
				message += "Could not find the directory of \"";
				message += _crystalProductionReportParameters.designedCrystalsDirectoryPath().generic_string();
				message += "\".";
			}

			throw System::IO::DirectoryNotFoundException{ typeid(*this), "setCrystalProductionReportParameters", message };
		}
	}


	_crystalDesignRecorder.setCrystalDesignRecordParameters(_crystalProductionReportParameters.crystalDesignRecordParameters());
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

void CrystalProductionReporter::outputFeasibleCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure) const
{
	OptimalCrystalStructure conventionalStructure = optimalCrystalStructure;
	conventionalStructure.conventionalizeStructure(_crystalProductionReportParameters.spaceGroupPrecision());
	{
		if (!(_crystalProductionReportParameters.needPiSymmetryCrystalData()) && (conventionalStructure.spaceGroupNumber() == ChemToolkit::Crystallography::Symmetry::SpaceGroupNumber{ 1 }))
			return;

		else
		{
			std::string structureFingerprint = conventionalStructure.toStructuralFingerprint();
			std::filesystem::path spaceGroupDirectoryPath = getSpaceGroupDirectoryPath(conventionalStructure.spaceGroupNumber(), optimalCrystalStructure.toChemicalComposition());

			constexpr size_type maxRepetitionCount = 50;
			std::lock_guard<std::mutex> guard{ s_reportMutex };


			for (size_type rep = 0; rep < maxRepetitionCount; ++rep)
			{
				try
				{
					auto stateAndDirectory = getOutputDirectoryPath(structureFingerprint, spaceGroupDirectoryPath);

					if (stateAndDirectory.first)
					{
						System::IO::Directory::createDirectory(stateAndDirectory.second, System::IO::Directory::CreateOptions::none);
						outputFeasibleCrystallographicData(conventionalStructure, optimalCrystalStructure, stateAndDirectory.second);


						std::filesystem::path fingerprintFilePath = stateAndDirectory.second;
						fingerprintFilePath /= s_fingerprintFilename;

						System::IO::FileStream fingerprintStreamWriter{ fingerprintFilePath, System::IO::FileStream::FileMode::createNew };
						fingerprintStreamWriter.write(structureFingerprint);
					}


					break;
				}

				catch (const System::ExceptionServices::IException&)
				{
					if ((1 + rep) < maxRepetitionCount)
						continue;
					else
						throw;
				}
			}
		}
	}
}

void CrystalProductionReporter::outputFeasibleCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure, const std::filesystem::path& producedDirectoryPath) const
{
	OptimalCrystalStructure conventionalStructure = optimalCrystalStructure;
	conventionalStructure.conventionalizeStructure(_crystalProductionReportParameters.spaceGroupPrecision());
	{
		if (!(_crystalProductionReportParameters.needPiSymmetryCrystalData()) && (conventionalStructure.spaceGroupNumber() == ChemToolkit::Crystallography::Symmetry::SpaceGroupNumber{ 1 }))
			return;

		else
		{
			std::string structureFingerprint = conventionalStructure.toStructuralFingerprint();
			std::filesystem::path spaceGroupDirectoryPath = getSpaceGroupDirectoryPath(conventionalStructure.spaceGroupNumber(), optimalCrystalStructure.toChemicalComposition());

			const size_type maxRepetitionCount = 50;
			std::lock_guard<std::mutex> guard{ s_reportMutex };


			for (size_type rep = 0; rep < maxRepetitionCount; ++rep)
			{
				try
				{
					auto stateAndDirectory = getOutputDirectoryPath(structureFingerprint, spaceGroupDirectoryPath);

					if (stateAndDirectory.first)
					{
						System::IO::Directory::createDirectory(stateAndDirectory.second, System::IO::Directory::CreateOptions::none);
						outputFeasibleCrystallographicData(conventionalStructure, optimalCrystalStructure, stateAndDirectory.second);


						std::filesystem::path fingerprintFilePath = stateAndDirectory.second;
						fingerprintFilePath /= s_fingerprintFilename;

						System::IO::FileStream fingerprintStreamWriter{ fingerprintFilePath, System::IO::FileStream::FileMode::createNew };
						fingerprintStreamWriter.write(structureFingerprint);


						std::filesystem::path productionReportDirectoryPath = stateAndDirectory.second;
						productionReportDirectoryPath /= getDirectoryName(producedDirectoryPath);
						System::IO::Directory::move(producedDirectoryPath, productionReportDirectoryPath);
					}

					else
					{
						std::filesystem::path productionReportDirectoryPath = stateAndDirectory.second;
						productionReportDirectoryPath /= getDirectoryName(producedDirectoryPath);
						System::IO::Directory::move(producedDirectoryPath, productionReportDirectoryPath);
					}


					break;
				}

				catch (const System::ExceptionServices::IException&)
				{
					if ((1 + rep) < maxRepetitionCount)
						continue;
					else
						throw;
				}
			}
		}
	}
}

void CrystalProductionReporter::outputInfeasibleCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure, const std::string& productionName) const
{
	if (_crystalProductionReportParameters.needInfeasibleCrystalData())
	{
		ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> chemicalComposition = optimalCrystalStructure.toChemicalComposition();

		std::filesystem::path outputDirectoryPath = _crystalProductionReportParameters.designedCrystalsDirectoryPath();
		outputDirectoryPath /= "Infeasible";
		outputDirectoryPath /= chemicalComposition.toChemicalSystem().toString();
		outputDirectoryPath /= chemicalComposition.toString();
		outputDirectoryPath /= productionName;


		std::lock_guard<std::mutex> guard{ s_reportMutex };
		System::IO::Directory::createDirectories(outputDirectoryPath, System::IO::Directory::CreateOptions::overwrite_existing);

		outputInfeasibleCrystallographicData(optimalCrystalStructure, outputDirectoryPath);
	}
}

void CrystalProductionReporter::outputInfeasibleCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure, const std::filesystem::path& producedDirectoryPath) const
{
	if (_crystalProductionReportParameters.needInfeasibleCrystalData())
	{
		ChemToolkit::Generic::ChemicalComposition<ChemToolkit::Generic::AtomicNumber> chemicalComposition = optimalCrystalStructure.toChemicalComposition();

		std::filesystem::path outputDirectoryPath = _crystalProductionReportParameters.designedCrystalsDirectoryPath();
		outputDirectoryPath /= "Infeasible";
		outputDirectoryPath /= chemicalComposition.toChemicalSystem().toString();
		outputDirectoryPath /= chemicalComposition.toString();


		std::lock_guard<std::mutex> guard{ s_reportMutex };
		System::IO::Directory::createDirectories(outputDirectoryPath, System::IO::Directory::CreateOptions::skip_existing);

		outputDirectoryPath /= getDirectoryName(producedDirectoryPath);
		{
			if (System::IO::Directory::exist(outputDirectoryPath))
				System::IO::Directory::deleteDirectory(outputDirectoryPath);
		}

		System::IO::Directory::move(producedDirectoryPath, outputDirectoryPath);

		outputInfeasibleCrystallographicData(optimalCrystalStructure, outputDirectoryPath);
	}
}

void CrystalProductionReporter::outputExceptionalCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure, const ChemicalComposition& chemicalComposition, const std::string& productionName) const
{
	if (_crystalProductionReportParameters.needExceptionalCrystalData())
	{
		std::filesystem::path outputDirectoryPath = _crystalProductionReportParameters.designedCrystalsDirectoryPath();
		outputDirectoryPath /= "Exception";
		outputDirectoryPath /= chemicalComposition.toChemicalSystem().toString();
		outputDirectoryPath /= chemicalComposition.toString();

		std::lock_guard<std::mutex> guard{ s_reportMutex };
		System::IO::Directory::createDirectories(outputDirectoryPath, System::IO::Directory::CreateOptions::skip_existing);

		outputDirectoryPath /= productionName;
		System::IO::Directory::createDirectory(outputDirectoryPath, System::IO::Directory::CreateOptions::overwrite_existing);


		if (!(optimalCrystalStructure.atoms().empty()))
			outputInfeasibleCrystallographicData(optimalCrystalStructure, outputDirectoryPath);
	}
}

void CrystalProductionReporter::outputExceptionalCrystalStructure(const OptimalCrystalStructure& optimalCrystalStructure, const ChemicalComposition& chemicalComposition, const std::filesystem::path& producedDirectoryPath) const
{
	if (_crystalProductionReportParameters.needExceptionalCrystalData())
	{
		std::filesystem::path outputDirectoryPath = _crystalProductionReportParameters.designedCrystalsDirectoryPath();
		outputDirectoryPath /= "Exception";
		outputDirectoryPath /= chemicalComposition.toChemicalSystem().toString();
		outputDirectoryPath /= chemicalComposition.toString();

		std::lock_guard<std::mutex> guard{ s_reportMutex };
		System::IO::Directory::createDirectories(outputDirectoryPath, System::IO::Directory::CreateOptions::skip_existing);

		outputDirectoryPath /= getDirectoryName(producedDirectoryPath);
		{
			if (System::IO::Directory::exist(outputDirectoryPath))
				System::IO::Directory::deleteDirectory(outputDirectoryPath);
		}


		System::IO::Directory::move(producedDirectoryPath, outputDirectoryPath);

		if (!(optimalCrystalStructure.atoms().empty()))
			outputInfeasibleCrystallographicData(optimalCrystalStructure, outputDirectoryPath);
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

std::filesystem::path CrystalProductionReporter::getSpaceGroupDirectoryPath(const SpaceGroupNumber spaceGroupNumber, const ChemicalComposition& composition) const
{
	std::filesystem::path spaceGroupDirectoryPath = _crystalProductionReportParameters.designedCrystalsDirectoryPath();
	{
		spaceGroupDirectoryPath /= "Optima";
		spaceGroupDirectoryPath /= composition.toChemicalSystem().toString();
		spaceGroupDirectoryPath /= composition.toString();
		spaceGroupDirectoryPath /= "SpaceGroup-";
		spaceGroupDirectoryPath += std::to_string(spaceGroupNumber);

		System::IO::Directory::createDirectories(spaceGroupDirectoryPath, System::IO::Directory::CreateOptions::skip_existing);
	}

	return spaceGroupDirectoryPath;
}

std::pair<bool, std::filesystem::path> CrystalProductionReporter::getOutputDirectoryPath(const std::string& outputFingerprint, const std::filesystem::path& spaceGroupPath) const
{
	size_type numDirectory = 1;
	{
		for (const auto& directoryPath : System::IO::Directory::enumerateDirectories(spaceGroupPath))
		{
			std::vector<std::filesystem::path> fingerprintFilePaths = System::IO::Directory::enumerateFiles(directoryPath, s_fingerprintFilename);


			if (fingerprintFilePaths.empty())
				throw System::IO::FileNotFoundException{ typeid(*this), "getStructureIndex", "Could not find fingerprint file." };

			else if (1 < fingerprintFilePaths.size())
				throw System::IO::InvalidFileException{ typeid(*this), "getStructureIndex", "There are multiple fingerprint files." };

			else
			{
				System::IO::FileStream fingerprintStreamReader{ fingerprintFilePaths.back(), System::IO::FileStream::FileMode::openRead };


				if (outputFingerprint == fingerprintStreamReader.readAllTexts())
					return std::make_pair(false, directoryPath);
				else
					++numDirectory;
			}
		}
	}


	std::filesystem::path outputDirectoryPath = spaceGroupPath;
	outputDirectoryPath /= "Type-";
	outputDirectoryPath += std::to_string(numDirectory);

	return std::make_pair(true, outputDirectoryPath);
}

void CrystalProductionReporter::outputFeasibleCrystallographicData(const OptimalCrystalStructure& conventionalCrystalStructure, const OptimalCrystalStructure& optimalCrystalStructure, const std::filesystem::path& outputDirectoryPath) const
{
	std::filesystem::path conventionalCifFilePath = outputDirectoryPath;
	std::filesystem::path optimalCifFilePath = outputDirectoryPath;
	{
		conventionalCifFilePath /= "conventionalStructure.cif";
		optimalCifFilePath /= "optimalStructure.cif";
	}

	ChemToolkit::Crystallography::IO::CifStreamWriter conventionalCifStreamWriter{ conventionalCifFilePath, System::IO::FileStream::FileMode::create };
	ChemToolkit::Crystallography::IO::CifStreamWriter optimalCifStreamWriter{ optimalCifFilePath, System::IO::FileStream::FileMode::create };

	conventionalCifStreamWriter.writeCrystallographicStructure(conventionalCrystalStructure);
	optimalCifStreamWriter.writeCrystallographicStructure(optimalCrystalStructure, 0.000001);
}

void CrystalProductionReporter::outputInfeasibleCrystallographicData(const OptimalCrystalStructure& optimalCrystalStructure, const std::filesystem::path& outputDirectoryPath) const
{
	std::filesystem::path optimalCifFilePath = outputDirectoryPath;
	optimalCifFilePath /= "optimalStructure.cif";

	ChemToolkit::Crystallography::IO::CifStreamWriter primitiveCifStreamWriter{ optimalCifFilePath, System::IO::FileStream::FileMode::create };
	primitiveCifStreamWriter.writeCrystallographicStructure(optimalCrystalStructure, 0.000001);
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
