#include "StreamReader.h"

#include <filesystem>
#include <sstream>
#include <regex>

#include "ArgumentOutOfRangeException.h"

#include <sstream>

using namespace System::IO;


// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// Constructors

StreamReader::StreamReader(std::vector<std::string>&& inputLines) noexcept
	: _allLines{ std::move(inputLines) }
{
}

StreamReader::StreamReader(const std::string& inputTexts)
	: _allLines{}
{
	std::istringstream inputStream{ inputTexts };

	for (std::string stringLine; std::getline(inputStream, stringLine);)
		_allLines.push_back(std::move(stringLine));
}

StreamReader::StreamReader(const StreamReader& sourceStreamReader, const std::string& startLine, const std::string& finalLine)
	: _allLines{}
{
	bool hasStartLine = false;
	{
		for (const auto& stringLine : sourceStreamReader._allLines)
		{
			if (hasStartLine)
			{
				_allLines.push_back(stringLine);

				if (stringLine == finalLine)
					break;
			}

			else
			{
				if (stringLine == startLine)
				{
					hasStartLine = true;
					_allLines.push_back(stringLine);
				}
			}
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
// Property

void StreamReader::setAllTexts(const std::string& inputTexts)
{
	_allLines.clear();
	std::istringstream inputStream{ inputTexts };

	for (std::string stringLine; std::getline(inputStream, stringLine);)
		_allLines.push_back(std::move(stringLine));
}

void StreamReader::pushBackAllTexts(const std::string& inputTexts)
{
	std::istringstream inputStream{ inputTexts };

	for (std::string stringLine; std::getline(inputStream, stringLine);)
		_allLines.push_back(std::move(stringLine));
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

StreamReader StreamReader::getListBlock(const std::string& prefixSymbol, const std::string& listName) const
{
	std::regex listNamePat;
	std::smatch listNameSm;
	{
		if (prefixSymbol == "&")
			listNamePat = std::regex{ "([\\_[:upper:]]+)" };
		else
			throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "getListBlock", "Prefix symbol is out of definitions." };
	}


	if (std::regex_match(listName, listNameSm, listNamePat))
	{
		std::regex spacePat{ "[\\s\\t]*" };
		std::string markerTexts;
		{
			markerTexts += prefixSymbol;
			markerTexts += listName;
		}

		StreamReader parameterStreamReader;
		{
			bool isFlagDiscovered = false;


			for (const auto& stringLine : _allLines)
			{
				if (isFlagDiscovered)
				{
					if (!(stringLine.empty()))
					{
						if (stringLine.front() == '&')
							break;

						else
						{
							std::smatch spaceSm;

							if (!(std::regex_match(stringLine, spaceSm, spacePat)))
								parameterStreamReader._allLines.push_back(stringLine);
						}
					}
				}

				else
				{
					if (stringLine == markerTexts)
						isFlagDiscovered = true;
				}
			}
		}

		return parameterStreamReader;
	}

	else
		return StreamReader{};
}

std::vector<std::pair<std::string, StreamReader>> StreamReader::enumerateListBlocks(const std::string& prefixSymbol) const
{
	std::vector<std::pair<std::string, StreamReader>> listBlocks;
	{
		std::regex listPat;
		std::smatch listSm;
		{
			if (prefixSymbol == "&")
				listPat = std::regex{ "[\\&]{1}([\\_[:upper:]]+)" };
			else
				throw System::ExceptionServices::ArgumentOutOfRangeException{ typeid(*this), "enumerateListBlocks", "Prefix symbol is out of definitions." };
		}


		for (const auto& stringLine : _allLines)
		{
			if (std::regex_match(stringLine, listSm, listPat))
				listBlocks.push_back(std::make_pair(stringLine, getListBlock(prefixSymbol, listSm.str(1))));
		}
	}

	return listBlocks;
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
// Static methods

int StreamReader::toIntegerValue(const std::string& inputText)
{
	std::regex pat("([\\-\\d]+)");
	std::smatch sm;

	if (std::regex_match(inputText, sm, pat))
		return std::stoi(sm.str(1));
	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::File::toIntegerValue", "Input text is invalid." };
}

double StreamReader::toFractionalValue(const std::string& inputText)
{
	std::regex integerPat{ "([\\-\\d]+)" };
	std::regex fractionPat{ "([\\-\\d]+)/([\\d]+)" };
	std::regex floatingPointPat{ "([e\\-\\.\\d]+)([\\(\\)\\d]+)*" };
	std::smatch sm;

	if (std::regex_match(inputText, sm, integerPat))
		return static_cast<double>(std::stoi(sm.str()));

	else if (std::regex_match(inputText, sm, fractionPat))
		return (static_cast<double>(std::stoi(sm.str(1))) / static_cast<double>(std::stoi(sm.str(2))));

	else if (std::regex_match(inputText, sm, floatingPointPat))
	{
		double value = 0.0;
		{
			std::stringstream valueStream{ sm.str(1) };
			valueStream >> value;
		}

		return value;
	}

	else
		throw System::ExceptionServices::ArgumentOutOfRangeException{ "System::IO::File::toFractionalValue", "Input text is invalid." };
}

std::vector<std::string> StreamReader::toParameterTuple(const std::string& stringLine)
{
	std::vector<std::string> parameterTuple;
	{
		std::regex pat("([\\:\\_\\-\\/\\+\\.[:alnum:]]+)");

		for (std::sregex_iterator iter(stringLine.begin(), stringLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			parameterTuple.push_back(iter->str(1));
	}

	return parameterTuple;
}

// Static methods
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

void StreamReader::pushBackStringLine(const std::string& inputLine)
{
	if (inputLine.find_first_of("\n") == std::string::npos)
		_allLines.push_back(inputLine);

	else
	{
		std::stringstream inputStream{ inputLine };


		for (std::string stringLine; std::getline(inputStream, stringLine);)
			_allLines.push_back(std::move(stringLine));
	}
}

void StreamReader::pushBackStringLine(std::string&& inputLine)
{
	if (inputLine.find_first_of("\n") == std::string::npos)
		_allLines.push_back(std::move(inputLine));

	else
	{
		std::stringstream inputStream{ inputLine };


		for (std::string stringLine; std::getline(inputStream, stringLine);)
			_allLines.push_back(std::move(stringLine));
	}
}

template <>
bool StreamReader::readValue<std::string>(const std::string& parameterLine, const std::string& dataFlag, std::string& outputText)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\-\\_\\/\\+\\.\\:[:alnum:]]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputText = sm.str(1);
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::filesystem::path>(const std::string& parameterLine, const std::string& dataFlag, std::filesystem::path& outputPath)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+[\"]{1}([\\s\\-\\_\\/\\+\\.\\:\\\\[:alnum:]]+){1}[\"]{1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputPath = std::filesystem::path{ sm.str(1) };
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::size_t>(const std::string& parameterLine, const std::string& dataFlag, std::size_t& outputInteger)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\d]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputInteger = static_cast<std::size_t>(std::stoi(sm.str(1)));
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<unsigned short>(const std::string& parameterLine, const std::string& dataFlag, unsigned short& outputInteger)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\d]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputInteger = static_cast<unsigned short>(std::stoi(sm.str(1)));
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<int>(const std::string& parameterLine, const std::string& dataFlag, int& outputInteger)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\+\\-\\d]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputInteger = std::stoi(sm.str(1));
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<short>(const std::string& parameterLine, const std::string& dataFlag, short& outputInteger)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\+\\-\\d]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		outputInteger = static_cast<short>(std::stoi(sm.str(1)));
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<double>(const std::string& parameterLine, const std::string& dataFlag, double& outputFloatingPointNumber)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([e\\-\\d\\.]+){1}([\\(\\)\\d]+)*[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		std::stringstream outputStream{ sm.str(1) };
		outputStream >> outputFloatingPointNumber;
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::pair<unsigned short, unsigned short>>(const std::string& parameterLine, const std::string& dataFlag, std::pair<unsigned short, unsigned short>& outputIntegerPair)
{
	std::string patString{ "[\\s\\t]*" };
	patString += dataFlag;
	patString += "[\\s\\t]+([\\d]+){1}[\\s\\t]+([\\d]+){1}[\\s]*";
	std::regex pat(patString);
	std::smatch sm;


	if (std::regex_match(parameterLine, sm, pat))
	{
		unsigned short firstValue = static_cast<unsigned short>(std::stoi(sm.str(1)));
		unsigned short secondValue = static_cast<unsigned short>(std::stoi(sm.str(2)));

		outputIntegerPair = std::make_pair(firstValue, secondValue);
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<std::string>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<std::string>& outputTexts)
{
	outputTexts.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\\:\\_\\-\\/\\+\\.[:alnum:]]+[\\s\\t]+)*([\\:\\_\\-\\/\\+\\.[:alnum:]]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([\\:\\_\\-\\/\\+\\.[:alnum:]]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputTexts.push_back(iter->str(1));

		outputTexts.erase(outputTexts.begin());
		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<std::filesystem::path>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<std::filesystem::path>& outputPaths)
{
	outputPaths.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\"]{1}[\\s\\:\\_\\-\\/\\+\\.\\\\[:alnum:]]+[\"]{1}[\\s\\t]+)*([\"]{1}[\\s\\:\\_\\-\\/\\+\\.\\\\[:alnum:]]+[\"]{1})[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("[\"]{1}([\\s\\:\\_\\-\\/\\+\\.\\\\[:alnum:]]+)[\"]{1}");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputPaths.push_back(std::filesystem::path{ iter->str(1) });

		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<std::size_t>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<std::size_t>& outputIntegers)
{
	outputIntegers.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\\d]+[\\s\\t]+)*([\\d]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([\\d]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputIntegers.push_back(static_cast<std::size_t>(std::stoi(iter->str(1))));

		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<unsigned short>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<unsigned short>& outputIntegers)
{
	outputIntegers.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\\d]+[\\s\\t]+)*([\\d]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([\\d]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputIntegers.push_back(static_cast<unsigned short>(std::stoi(iter->str(1))));

		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<short>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<short>& outputIntegers)
{
	outputIntegers.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\\+\\-\\d]+[\\s\\t]+)*([\\+\\-\\d]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([\\-\\d]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputIntegers.push_back(static_cast<short>(std::stoi(iter->str(1))));

		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<int>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<int>& outputIntegers)
{
	outputIntegers.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([\\+\\-\\d]+[\\s\\t]+)*([\\+\\-\\d]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([\\-\\d]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
			outputIntegers.push_back(std::stoi(iter->str(1)));

		return true;
	}

	else
		return false;
}

template <>
bool StreamReader::readValue<std::vector<double>>(const std::string& parameterLine, const std::string& dataFlag, std::vector<double>& outputFloatingPointNumbers)
{
	outputFloatingPointNumbers.clear();

	std::string totalPatString{ "[\\s\\t]*" };
	totalPatString += dataFlag;
	totalPatString += "[\\s\\t]+([e\\-\\d\\.]+[\\s\\t]+)*([e\\-\\d\\.]+)[\\s]*";
	std::regex totalPat(totalPatString);
	std::smatch totalSm;


	if (std::regex_match(parameterLine, totalSm, totalPat))
	{
		std::regex pat("([e\\-\\d\\.]+)");

		for (std::sregex_iterator iter(parameterLine.begin(), parameterLine.end(), pat); iter != std::sregex_iterator{}; ++iter)
		{
			double value = 0.0;
			{
				std::stringstream outputStream{ iter->str(1) };
				outputStream >> value;
			}

			outputFloatingPointNumbers.push_back(value);
		}

		return true;
	}

	else
		return false;
}

// Private methods
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
// **********************************************************************************************************************************************************************************************************************************************************************************************
