#pragma once

struct Orf;
void OrfScanning(std::ostream& os, const std::string& seq, int len);
bool IsInitiation(std::string coden);
bool IsTermination(std::string coden);

int main();
