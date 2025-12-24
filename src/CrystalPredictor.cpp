#include "CrystalPredictor.h"

#include <mpi.h>
#include <iostream>
#include <vector>
#include <thread>

#include "ThreadingPolicy.h"
#include "MpiPolicy.h"
#include "DateTime.h"

#include "Directory.h"
#include "DirectoryNotFoundException.h"
#include "File.h"
#include "FileNotFoundException.h"
#include "LogStreamWriter.h"
#include "StreamWriter.h"

#include "ConstrainingAtom.h"
#include "ConstrainingCrystalStructure.h"
#include "OptimalCrystalStructure.h"

#include "FeasiblePolyhedraConnectionsDictionary.h"

using namespace MathematicalCrystalChemistry::Prediction;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

CrystalPredictor::size_type CrystalPredictor::ProduceCrystals::s_structureDesigning;
CrystalPredictor::size_type CrystalPredictor::ProduceCrystals::s_maxStructureDesigning;
CrystalPredictor::ChemicalComposition CrystalPredictor::ProduceCrystals::s_chemicalComposition;
std::mutex CrystalPredictor::ProduceCrystals::s_productionMutex;


CrystalPredictor::ProduceCrystals::ProduceCrystals(const size_type rank) noexcept
	: _stdFilePath{}
	, _threadRank{ rank }
	, _randomStructureGenerator{}
	, _crystalDesigner{}
	, _crystalProductionReporter{}
	, m_sampleIndex{ 0 }
{
}

CrystalPredictor::CrystalPredictor() noexcept
	: _isStdoutEnabled{ true }
	, _crystalPredictionTask{}
	, m_stdFilePath{}
	, m_crystalProducers{}
{
}

CrystalPredictor::CrystalPredictor(const CrystalPredictionTask& task)
	: _isStdoutEnabled{ true }
	, _crystalPredictionTask{ task }
	, m_stdFilePath{}
	, m_crystalProducers{}
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

void CrystalPredictor::ProduceCrystals::operator()()
{
	try
	{
		System::IO::LogStreamWriter logStreamWriter{ _crystalProductionReporter.crystalDesignRecorder().mpiCrystalProductionDirectoryPath() };

		MathematicalCrystalChemistry::CrystalModel::Components::ConstrainingCrystalStructure constrainingCrystalStructure;
		_randomStructureGenerator.setGeneratingChemicalComposition(s_chemicalComposition);

		const size_type numDesigningUnit = (s_maxStructureDesigning / 10);


		while (shouldDesign())
		{
			{
				if ((s_structureDesigning % numDesigningUnit) == 0)
					reportCeaselessGeneration();
			}
			std::filesystem::path producedDirectoryPath = _crystalProductionReporter.crystalDesignRecorder().mpiCrystalProductionDirectoryPath();
			const std::string productionName = getProductionName();
			producedDirectoryPath /= productionName;

			MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure optimalCrystalStructure;


			try
			{
				_randomStructureGenerator.next(constrainingCrystalStructure);
				{
					if (_crystalProductionReporter.crystalDesignRecorder().needMdRecord())
					{
						_crystalProductionReporter.crystalDesignRecorder().setProductionName(productionName);
						_crystalDesigner.execute(constrainingCrystalStructure, _crystalProductionReporter.crystalDesignRecorder());
					}

					else
						_crystalDesigner.execute(constrainingCrystalStructure);
				}
				constrainingCrystalStructure.setFeasibleErrorRate(_crystalDesigner.preciseStructuralOptimizationParameters().feasibleGeometricalConstraintErrorRate());


				optimalCrystalStructure = MathematicalCrystalChemistry::CrystalModel::Components::OptimalCrystalStructure{ constrainingCrystalStructure };
				{
					if (_crystalProductionReporter.crystalDesignRecorder().needMdRecord())
					{
						if (constrainingCrystalStructure.isFeasible())
							_crystalProductionReporter.outputFeasibleCrystalStructure(optimalCrystalStructure, producedDirectoryPath);
						else
							_crystalProductionReporter.outputInfeasibleCrystalStructure(optimalCrystalStructure, producedDirectoryPath);
					}

					else
					{
						if (constrainingCrystalStructure.isFeasible())
							_crystalProductionReporter.outputFeasibleCrystalStructure(optimalCrystalStructure);
						else
							_crystalProductionReporter.outputInfeasibleCrystalStructure(optimalCrystalStructure, productionName);
					}
				}
			}


			catch (const System::ExceptionServices::IException& e)
			{
				{
					if (_stdFilePath.empty())
						logStreamWriter.write(e);
				}

				if (_crystalProductionReporter.crystalDesignRecorder().needMdRecord())
					_crystalProductionReporter.outputExceptionalCrystalStructure(optimalCrystalStructure, toBasicChemicalComposition(s_chemicalComposition), producedDirectoryPath);
				else
					_crystalProductionReporter.outputExceptionalCrystalStructure(optimalCrystalStructure, toBasicChemicalComposition(s_chemicalComposition), productionName);
			}

			catch (const std::exception& e)
			{
				{
					if (_stdFilePath.empty())
						logStreamWriter.write(e);
				}

				if (_crystalProductionReporter.crystalDesignRecorder().needMdRecord())
					_crystalProductionReporter.outputExceptionalCrystalStructure(optimalCrystalStructure, toBasicChemicalComposition(s_chemicalComposition), producedDirectoryPath);
				else
					_crystalProductionReporter.outputExceptionalCrystalStructure(optimalCrystalStructure, toBasicChemicalComposition(s_chemicalComposition), productionName);
			}
		}

		_crystalDesigner.setProductChemicalComposition();
	}


	catch (const System::ExceptionServices::IException& e)
	{
		if (_stdFilePath.empty())
			std::cout << e.toString() << std::endl;

		/*
		else
		{
			System::IO::FileStream stdFileStream{ _stdFilePath, System::IO::FileStream::FileMode::append };
			stdFileStream.write(e.toString());
		}
		*/
	}

	catch (const std::exception& e)
	{
		if (_stdFilePath.empty())
			std::cout << e.what() << std::endl;

		/*
		else
		{
			System::IO::FileStream stdFileStream{ _stdFilePath, System::IO::FileStream::FileMode::append };
			stdFileStream.write(e.what());
		}
		*/
	}
}

void CrystalPredictor::execute() const
{
	std::filesystem::path crystalProductionDirectoryPath{ _crystalPredictionTask.crystalProductionReportParameters().crystalDesignRecordParameters().getCrystalProductionDirectoryPath() };
	MPI_Barrier(MPI_COMM_WORLD);

	createMpiProductionDirectoryPaths(crystalProductionDirectoryPath);
	updateStdFilePath(crystalProductionDirectoryPath);

	try
	{
		createLogFiles(crystalProductionDirectoryPath);
		MPI_Barrier(MPI_COMM_WORLD);

		initializeCrystalProducers(crystalProductionDirectoryPath);
		validateProductionPaths(crystalProductionDirectoryPath);
		MPI_Barrier(MPI_COMM_WORLD);


		for (const auto& compositionAndGenerating : _crystalPredictionTask.crystalDesignRequest())
		{
			ProduceCrystals::initializeStructureProducing();
			ProduceCrystals::setMaxStructureProducing(getMaxCrystalProducing(compositionAndGenerating.second));
			ProduceCrystals::setChemicalComposition(compositionAndGenerating.first);

			MathematicalCrystalChemistry::CrystalModel::Constraints::FeasiblePolyhedraConnectionsDictionary::initialize(compositionAndGenerating.first);


			std::vector<std::thread> operatingSystems;
			{
				for (std::size_t threadRank = 0; threadRank < System::Parallel::ThreadingPolicy::maxThreading(); ++threadRank)
					operatingSystems.push_back(std::thread{ m_crystalProducers[threadRank] });
			}

			for (auto& operatingSystem : operatingSystems)
				operatingSystem.join();


			reportJobFinalization(compositionAndGenerating.first, getMaxCrystalProducing(compositionAndGenerating.second));
		}


		reportAllJobFinalization();
	}

	catch (const System::ExceptionServices::IException& e)
	{
		if (_isStdoutEnabled)
			std::cout << e.toString() << std::endl;

		/*
		else
		{
			System::IO::FileStream stdFileStream{ m_stdFilePath, System::IO::FileStream::FileMode::append };
			stdFileStream.write(e.toString());
		}
		*/
	}

	catch (const std::exception& e)
	{
		if (_isStdoutEnabled)
			std::cout << e.what() << std::endl;

		/*
		else
		{
			System::IO::FileStream stdFileStream{ m_stdFilePath, System::IO::FileStream::FileMode::append };
			stdFileStream.write(e.what());
		}
		*/
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

void CrystalPredictor::ProduceCrystals::reportCeaselessGeneration() const
{
	System::IO::StreamWriter streamWriter;
	{
		double designedPercentage = 100.0;
		designedPercentage *= (static_cast<double>(s_structureDesigning) / static_cast<double>(s_maxStructureDesigning));

		streamWriter.write(System::Utility::DateTime::toTimeString(System::Utility::DateTime::getNowTime()));
		streamWriter.write(": Process_");
		streamWriter.write(System::Parallel::MpiPolicy::mpiRank());
		streamWriter.write(" has designed ");
		streamWriter.write(designedPercentage, 2);
		streamWriter.write("% structures.");
		streamWriter.breakLine();
	}


	if (_stdFilePath.empty())
		std::cout << streamWriter.allTexts();

	/*
	else
	{
		System::IO::FileStream stdFileStream{ _stdFilePath, System::IO::FileStream::FileMode::append };
		stdFileStream.write(streamWriter.allTexts());
	}
	*/
}

std::string CrystalPredictor::ProduceCrystals::getProductionName() const
{
	std::string name;
	{
		name += "Mpi_";
		name += std::to_string(System::Parallel::MpiPolicy::mpiRank());
		name += "_Production_";
		name += std::to_string(m_sampleIndex);
	}

	return name;
}

void CrystalPredictor::updateStdFilePath(const std::filesystem::path& crystalProductionDirectoryPath) const
{
	if (!_isStdoutEnabled)
	{
		m_stdFilePath = crystalProductionDirectoryPath;
		m_stdFilePath /= System::Parallel::MpiPolicy::getMpiDirectoryName();
		m_stdFilePath /= "stdout.txt";
	}
}

void CrystalPredictor::createMpiProductionDirectoryPaths(const std::filesystem::path& crystalProductionDirectoryPath) const
{
	if (System::Parallel::MpiPolicy::mpiRank() == 0)
	{
		System::IO::Directory::createDirectory(crystalProductionDirectoryPath, System::IO::Directory::CreateOptions::overwrite_existing);


		for (size_type rank = 0; rank < System::Parallel::MpiPolicy::mpiProcessing(); ++rank)
		{
			std::filesystem::path mpiCrystalProductionDirectoryPath{ crystalProductionDirectoryPath };
			mpiCrystalProductionDirectoryPath /= "Mpi_";
			mpiCrystalProductionDirectoryPath += std::to_string(rank);

			System::IO::Directory::createDirectory(mpiCrystalProductionDirectoryPath, System::IO::Directory::CreateOptions::none);
		}
	}
}

void CrystalPredictor::createLogFiles(const std::filesystem::path& crystalProductionDirectoryPath) const
{
	if (_isStdoutEnabled)
	{
		if (System::Parallel::MpiPolicy::mpiRank() == 0)
		{
			for (size_type rank = 0; rank < System::Parallel::MpiPolicy::mpiProcessing(); ++rank)
			{
				std::filesystem::path mpiCrystalProductionDirectoryPath{ crystalProductionDirectoryPath };
				mpiCrystalProductionDirectoryPath /= "Mpi_";
				mpiCrystalProductionDirectoryPath += std::to_string(rank);

				System::IO::LogStreamWriter logStreamWriter{ mpiCrystalProductionDirectoryPath };
			}
		}
	}
}

void CrystalPredictor::initializeCrystalProducers(const std::filesystem::path& crystalProductionDirectoryPath) const
{
	std::filesystem::path mpiCrystalProductionDirectoryPath{ crystalProductionDirectoryPath };
	mpiCrystalProductionDirectoryPath /= "Mpi_";
	mpiCrystalProductionDirectoryPath += std::to_string(System::Parallel::MpiPolicy::mpiRank());

	m_crystalProducers.clear();
	{
		if (_isStdoutEnabled)
		{
			for (std::size_t threadRank = 0; threadRank < System::Parallel::ThreadingPolicy::maxThreading(); ++threadRank)
			{
				ProduceCrystals produceCrystals{ threadRank };
				produceCrystals.setCrystalDesignParameters(_crystalPredictionTask.crystalDesignParameters());
				produceCrystals.setCrystalProductionReportParameters(_crystalPredictionTask.crystalProductionReportParameters(), mpiCrystalProductionDirectoryPath);

				m_crystalProducers.push_back(std::move(produceCrystals));
			}
		}

		else
		{
			for (std::size_t threadRank = 0; threadRank < System::Parallel::ThreadingPolicy::maxThreading(); ++threadRank)
			{
				ProduceCrystals produceCrystals{ threadRank };
				produceCrystals.setStdFilePath(m_stdFilePath);
				produceCrystals.setCrystalDesignParameters(_crystalPredictionTask.crystalDesignParameters());
				produceCrystals.setCrystalProductionReportParameters(_crystalPredictionTask.crystalProductionReportParameters(), mpiCrystalProductionDirectoryPath);

				m_crystalProducers.push_back(std::move(produceCrystals));
			}
		}
	}
}

void CrystalPredictor::validateProductionPaths(const std::filesystem::path& crystalProductionDirectoryPath) const
{
	std::filesystem::path mpiCrystalProductionDirectoryPath{ crystalProductionDirectoryPath };
	mpiCrystalProductionDirectoryPath /= "Mpi_";
	mpiCrystalProductionDirectoryPath += std::to_string(System::Parallel::MpiPolicy::mpiRank());

	std::filesystem::path mpiLogFilePath{ mpiCrystalProductionDirectoryPath };
	mpiLogFilePath /= System::IO::LogStreamWriter::getLogFilename();


	if (!(System::IO::Directory::exist(crystalProductionDirectoryPath)))
	{
		std::string errorMessage;
		{
			errorMessage += "Could not find the directory of \"";
			errorMessage += crystalProductionDirectoryPath.generic_string();
			errorMessage += "\".";
		}

		throw System::IO::DirectoryNotFoundException{ typeid(*this), "validateProductionPaths", errorMessage };
	}

	if (!(System::IO::Directory::exist(mpiCrystalProductionDirectoryPath)))
	{
		std::string errorMessage;
		{
			errorMessage += "Could not find the directory of \"";
			errorMessage += mpiCrystalProductionDirectoryPath.generic_string();
			errorMessage += "\".";
		}

		throw System::IO::DirectoryNotFoundException{ typeid(*this), "validateProductionPaths", errorMessage };
	}

	if (_isStdoutEnabled && !(System::IO::File::exist(mpiLogFilePath)))
	{
		std::string errorMessage;
		{
			errorMessage += "Could not find the file of \"";
			errorMessage += mpiLogFilePath.generic_string();
			errorMessage += "\".";
		}

		throw System::IO::FileNotFoundException{ typeid(*this), "validateProductionPaths", errorMessage };
	}
}

void CrystalPredictor::reportAllJobFinalization() const
{
	std::string message;
	{
		message += "\n================================================================================================================\n";;
		message += "================================================================================================================\n";;
		message += "\tFinish all the creation of crystal structure prototypes.\n\n\n";
	}


	if (_isStdoutEnabled)
		std::cout << message << std::endl;
	else
	{
		System::IO::FileStream stdFileStream{ m_stdFilePath, System::IO::FileStream::FileMode::append };
		stdFileStream.write(message);
	}
}

void CrystalPredictor::reportJobFinalization(const ChemicalComposition& chemicalComposition, const size_type structureGenerating) const
{
	std::string message;
	{
		message += "MPI ";
		message += std::to_string(System::Parallel::MpiPolicy::mpiRank());
		message += ":  Finish ";
		message += std::to_string(structureGenerating);
		message += " generating of ";
		message += chemicalComposition.toString();
	}


	if (_isStdoutEnabled)
		std::cout << message << std::endl;

	/*
	else
	{
		System::IO::FileStream stdFileStream{ m_stdFilePath, System::IO::FileStream::FileMode::append };
		stdFileStream.write(message);
	}
	*/
}

CrystalPredictor::size_type CrystalPredictor::getMaxCrystalProducing(const size_type requestedCrystalProducing) const
{
	size_type structureProducing = (requestedCrystalProducing / System::Parallel::MpiPolicy::mpiProcessing());
	{
		size_type remainder = (requestedCrystalProducing % System::Parallel::MpiPolicy::mpiProcessing());

		if (System::Parallel::MpiPolicy::mpiRank() < remainder)
			structureProducing += 1;
	}

	return structureProducing;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
