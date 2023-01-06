#pragma once

#include <random>

class RandomNumber
{
public:
    static double Uniform()
    {
      return getInstance().uniformDistribution_(getInstance().mt);
    }
    static double Gaussian()
    {
      return getInstance().normalDistribution_(getInstance().mt);
    }
    static size_t Integer(size_t i, size_t j)
    {
      return i + static_cast<size_t>(static_cast<double>(j + 1 - i) * Uniform());
    }
    static uint64_t UInt64()
    {
      return getInstance().uniformUInt64Distribution_(getInstance().mt);
    }
private:
    RandomNumber()
    {
      std::random_device rd;
      mt = std::mt19937_64(rd());
      //mt = std::mt19937_64(10);
      uniformDistribution_ = std::uniform_real_distribution<double>(0.0, 1.0);
      uniformUInt64Distribution_ = std::uniform_int_distribution<uint64_t>();
      normalDistribution_ = std::normal_distribution<double>();
    }
    ~RandomNumber() {}
    
    static RandomNumber& getInstance() {
        static RandomNumber s;
        return s;
    }

    RandomNumber(RandomNumber const&) = delete;
    RandomNumber& operator= (RandomNumber const&) = delete;

    std::mt19937_64 mt;
    std::uniform_real_distribution<double> uniformDistribution_;
    std::normal_distribution<double> normalDistribution_;
    std::uniform_int_distribution<uint64_t> uniformUInt64Distribution_;
};
