#include "Directory.h"

#include <regex>

using namespace System::IO;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Methods

void Directory::copy(const std::filesystem::path& sourceDirectoryPath, const std::filesystem::path& destinationDirectoryPath, const CopyOptions& copyOption)
{
	validateCreatingDirectoryPath(destinationDirectoryPath);


	switch (copyOption)
	{
	case CopyOptions::none:
		std::filesystem::copy(sourceDirectoryPath, destinationDirectoryPath, std::filesystem::copy_options::none);
		break;

	case CopyOptions::skip_existing:
		std::filesystem::copy(sourceDirectoryPath, destinationDirectoryPath, std::filesystem::copy_options::skip_existing);
		break;

	case CopyOptions::overwrite_existing:
		std::filesystem::copy(sourceDirectoryPath, destinationDirectoryPath, std::filesystem::copy_options::overwrite_existing);
		break;

	case CopyOptions::recursive:
		std::filesystem::copy(sourceDirectoryPath, destinationDirectoryPath, std::filesystem::copy_options::recursive);
		break;

	case CopyOptions::directories_only:
		std::filesystem::copy(sourceDirectoryPath, destinationDirectoryPath, std::filesystem::copy_options::directories_only);
		break;

	default:
		throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::copy", "\"copyOption\" is invalid." };
	}
}

void Directory::move(const std::filesystem::path& sourceDirectoryPath, const std::filesystem::path& destinationDirectoryPath, const CreateOptions& createOption)
{
	validateCreatingDirectoryPath(destinationDirectoryPath);


	if (std::filesystem::exists(destinationDirectoryPath))
	{
		switch (createOption)
		{
		case CreateOptions::none:
			throw IOException{ "System::IO::Directory::move", "Destination directory path is the path to an existing object." };

		case CreateOptions::skip_existing:
			return;

		case CreateOptions::overwrite_existing:
			std::filesystem::remove_all(destinationDirectoryPath);
			break;

		default:
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::move", "\"createOption\" is invalid." };
		}
	}


	std::filesystem::rename(sourceDirectoryPath, destinationDirectoryPath);
}

void Directory::createDirectory(const std::filesystem::path& directoryPath, const CreateOptions& createOption)
{
	validateCreatingDirectoryPath(directoryPath);


	if (std::filesystem::exists(directoryPath))
	{
		switch (createOption)
		{
		case CreateOptions::none:
		{
			std::string message;
			{
				message += "Directory path of \"";
				message += directoryPath.generic_string();
				message += "\" is the path to an existing object.";
			}

			throw IOException{ "System::IO::Directory::createDirectory", message };
		}

		case CreateOptions::skip_existing:
			return;

		case CreateOptions::overwrite_existing:
			std::filesystem::remove_all(directoryPath);
			break;

		default:
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::createDirectory", "\"createOption\" is invalid." };
		}
	}


	if (!(std::filesystem::create_directory(directoryPath)))
	{
		std::string errorMessage;
		{
			errorMessage += "Could not create the directory of \"";
			errorMessage += directoryPath.generic_string();
			errorMessage += "\".";
		}

		throw IOException{ "System::IO::Directory::createDirectory", errorMessage };
	}
}

void Directory::createDirectories(const std::filesystem::path& directoryPath, const CreateOptions& createOption)
{
	if (std::filesystem::exists(directoryPath))
	{
		switch (createOption)
		{
		case CreateOptions::none:
			throw IOException{ "System::IO::Directory::createDirectories", "Directory path is the path to an existing object." };

		case CreateOptions::skip_existing:
			return;

		case CreateOptions::overwrite_existing:
			std::filesystem::remove_all(directoryPath);
			break;

		default:
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::createDirectories", "\"createOption\" is invalid." };
		}
	}


	if (!(std::filesystem::create_directories(directoryPath)))
	{
		std::string errorMessage;
		{
			errorMessage += "Could not create the directories of \"";
			errorMessage += directoryPath.generic_string();
			errorMessage += "\".";
		}

		throw System::IO::IOException{ "System::IO::Directory::createDirectories", errorMessage };
	}
}

std::size_t System::IO::Directory::countDirectories(const std::filesystem::path& directoryPath, const std::string& directoryNamePattern, const SearchOptions& searchOption)
{
	std::size_t numDirectories = 0;
	{
		std::regex pat(directoryNamePattern);
		std::smatch sm;


		if (searchOption == SearchOptions::TopDirectoryOnly)
		{
			for (const std::filesystem::directory_entry& directoryEntry : std::filesystem::directory_iterator(directoryPath))
			{
				if (directoryEntry.is_directory())
				{
					std::string directoryName = getDirectoryName(directoryEntry.path());

					if (std::regex_match(directoryName, sm, pat))
						++numDirectories;
				}
			}
		}

		else if (searchOption == SearchOptions::AllDirectories)
		{
			for (const std::filesystem::directory_entry& directoryEntry : std::filesystem::recursive_directory_iterator(directoryPath))
			{
				if (directoryEntry.is_directory())
				{
					std::string directoryName = getDirectoryName(directoryEntry.path());

					if (std::regex_match(directoryName, sm, pat))
						++numDirectories;
				}
			}
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::enumerateDirectories", "\"searchOption\" is invalid." };
	}

	return numDirectories;
}

std::vector<std::filesystem::path> Directory::enumerateDirectories(const std::filesystem::path& directoryPath, const std::string& directoryNamePattern, const SearchOptions& searchOption)
{
	std::vector<std::filesystem::path> directoryPaths;
	{
		std::regex pat(directoryNamePattern);
		std::smatch sm;


		if (searchOption == SearchOptions::TopDirectoryOnly)
		{
			for (const std::filesystem::directory_entry& directoryEntry : std::filesystem::directory_iterator(directoryPath))
			{
				if (directoryEntry.is_directory())
				{
					std::string directoryName = getDirectoryName(directoryEntry.path());

					if (std::regex_match(directoryName, sm, pat))
						directoryPaths.push_back(directoryEntry.path());
				}
			}
		}

		else if (searchOption == SearchOptions::AllDirectories)
		{
			for (const std::filesystem::directory_entry& directoryEntry : std::filesystem::recursive_directory_iterator(directoryPath))
			{
				if (directoryEntry.is_directory())
				{
					std::string directoryName = getDirectoryName(directoryEntry.path());

					if (std::regex_match(directoryName, sm, pat))
						directoryPaths.push_back(directoryEntry.path());
				}
			}
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::enumerateDirectories", "\"searchOption\" is invalid." };
	}

	return directoryPaths;
}

std::vector<std::filesystem::path> Directory::enumerateFiles(const std::filesystem::path& directoryPath, const std::string& fileNamePattern, const SearchOptions& searchOption)
{
	std::vector<std::filesystem::path> filePaths;
	{
		std::regex pat(fileNamePattern);
		std::smatch sm;


		if (searchOption == SearchOptions::TopDirectoryOnly)
		{
			for (const std::filesystem::directory_entry& fileEntry : std::filesystem::directory_iterator(directoryPath))
			{
				if (fileEntry.is_regular_file())
				{
					std::string fileName = fileEntry.path().filename().string();

					if (std::regex_match(fileName, sm, pat))
						filePaths.push_back(fileEntry.path());
				}
			}
		}

		else if (searchOption == SearchOptions::AllDirectories)
		{
			for (const std::filesystem::directory_entry& fileEntry : std::filesystem::recursive_directory_iterator(directoryPath))
			{
				if (fileEntry.is_regular_file())
				{
					std::string fileName = fileEntry.path().filename().string();

					if (std::regex_match(fileName, sm, pat))
						filePaths.push_back(fileEntry.path());
				}
			}
		}

		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::Directory::enumerateDirectories", "\"searchOption\" is invalid." };
	}

	return filePaths;
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

void Directory::validateDirectory(const std::filesystem::path& directoryPath)
{
	if (!(exist(directoryPath)))
	{
		std::string errorMessage;
		{
			errorMessage += "\"";
			errorMessage += directoryPath.generic_string();
			errorMessage += "\" is invalid.";
		}

		throw FileNotFoundException{ "System::IO::Directory::validateDirectory", errorMessage };
	}
}

void Directory::validateCreatingDirectoryPath(const std::filesystem::path& directoryPath)
{
	if (!(std::filesystem::is_directory(directoryPath.parent_path())))
	{
		std::string errorMessage;
		{
			errorMessage += "Parent path of \"";
			errorMessage += directoryPath.generic_string();
			errorMessage += "\" is not a path to a directory.";
		}

		throw DirectoryNotFoundException{ "System::IO::Directory::validateCreatingDirectoryPath", errorMessage };
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
