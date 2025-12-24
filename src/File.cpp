#include "File.h"

#include <fstream>

#include "FileNotFoundException.h"

#include "StreamReader.h"
#include "StreamWriter.h"

using namespace System::IO;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Methods

void File::copy(const std::filesystem::path& sourceFilePath, const std::filesystem::path& destinationFilePath, const CreateOptions& createOption)
{
	validateCreatingFilePath(destinationFilePath);


	switch (createOption)
	{
	case CreateOptions::none:
		std::filesystem::copy_file(sourceFilePath, destinationFilePath, std::filesystem::copy_options::none);
		break;

	case CreateOptions::skip_existing:
		std::filesystem::copy_file(sourceFilePath, destinationFilePath, std::filesystem::copy_options::skip_existing);
		break;

	case CreateOptions::overwrite_existing:
		std::filesystem::copy_file(sourceFilePath, destinationFilePath, std::filesystem::copy_options::overwrite_existing);
		break;

	default:
		throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::File::copy", "\"createOption\" is invalid." };
	}
}

void File::copyTo(const std::filesystem::path& sourceFilePath, const std::filesystem::path& destinationDirectoryPath, const CreateOptions& createOption)
{
	if (std::filesystem::is_directory(destinationDirectoryPath))
	{
		std::filesystem::path destinationFilePath = destinationDirectoryPath;
		destinationFilePath /= sourceFilePath.filename();
		copy(sourceFilePath, destinationFilePath, createOption);
	}

	else
		throw DirectoryNotFoundException{ "System::IO::File::copyTo", "Destination directory path is not a path to a directory." };
}

void File::move(const std::filesystem::path& sourceFilePath, const std::filesystem::path& destinationFilePath, const CreateOptions& createOption)
{
	validateCreatingFilePath(destinationFilePath);


	if (std::filesystem::exists(destinationFilePath))
	{
		switch (createOption)
		{
		case CreateOptions::none:
			throw IOException{ "System::IO::File::move", "Destination file path is the path to an existing object." };

		case CreateOptions::skip_existing:
			return;

		case CreateOptions::overwrite_existing:
			std::filesystem::remove_all(destinationFilePath);
			break;

		default:
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::File::move", "\"createOption\" is invalid." };
		}
	}

	std::filesystem::rename(sourceFilePath, destinationFilePath);
}

void File::moveTo(const std::filesystem::path& sourceFilePath, const std::filesystem::path& destinationDirectoryPath, const CreateOptions& createOption)
{
	if (std::filesystem::is_directory(destinationDirectoryPath))
	{
		std::filesystem::path destinationFilePath = destinationDirectoryPath;
		destinationFilePath /= sourceFilePath.filename();
		move(sourceFilePath, destinationFilePath, createOption);
	}

	else
		throw DirectoryNotFoundException{ "System::IO::File::moveTo", "Destination directory path is not a path to a directory." };
}

void File::createFile(const std::filesystem::path& filePath, const CreateOptions& createOption)
{
	validateCreatingFilePath(filePath);

	if (std::filesystem::exists(filePath))
	{
		switch (createOption)
		{
		case CreateOptions::none:
			throw System::IO::IOException{ "System::IO::File::createFile", "File path is the path to an existing object." };

		case CreateOptions::skip_existing:
			return;

		case CreateOptions::overwrite_existing:
			std::filesystem::remove_all(filePath);
			break;

		default:
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::File::createFile", "\"createOption\" is invalid." };
		}
	}


	std::ofstream outputFileStream;
	outputFileStream.open(filePath, std::ios_base::out);
	{
		if (outputFileStream)
			outputFileStream.close();
		else
		{
			std::string errorMessage;
			{
				errorMessage += "Could not open the file of \"";
				errorMessage += filePath.generic_string();
				errorMessage += "\".";
			}

			throw System::IO::FileNotFoundException{ "System::IO::File::createFile", errorMessage };
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

void File::validateRegularFile(const std::filesystem::path& filePath)
{
	if (!(exist(filePath)))
	{
		std::string errorMessage;
		{
			errorMessage += "\"";
			errorMessage += filePath.generic_string();
			errorMessage += "\" is invalid.";
		}

		throw FileNotFoundException{ "System::IO::File::validateRegularFile", errorMessage };
	}
}

void File::validateCreatingFilePath(const std::filesystem::path& filePath)
{
	if (filePath.has_filename())
	{
		if (!(std::filesystem::is_directory(filePath.parent_path())))
		{
			std::string errorMessage;
			{
				errorMessage += "Parent path of \"";
				errorMessage += filePath.parent_path().generic_string();
				errorMessage += "\" is not a path to a directory.";
			}

			throw DirectoryNotFoundException{ "System::IO::File::validateCreatingFilePath", errorMessage };
		}
	}

	else
	{
		std::string errorMessage;
		{
			errorMessage += "\"";
			errorMessage += filePath.generic_string();
			errorMessage += "\" does not constain filename.";
		}

		throw IOException{ "System::IO::File::validateCreatingFilePath", errorMessage };
	}
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
