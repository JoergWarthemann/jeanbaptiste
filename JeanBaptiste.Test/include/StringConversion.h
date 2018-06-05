#pragma once

#include <codecvt>
#include <locale>
#include <string>

namespace Utilities {

static std::string narrow(const std::wstring& wideString)
{
    std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
    return converter.to_bytes(wideString);
}

static std::wstring widen(const std::string& narrowString)
{
    std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
    return converter.from_bytes(narrowString);
}

}
