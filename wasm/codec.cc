#include <string>

#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>

using namespace brunsli;

extern "C" {

std::string* BrunsliToJpeg(const uint8_t* data, size_t length) {
  std::string* result = nullptr;
  JPEGData jpg;
  BrunsliStatus status =
      BrunsliDecodeJpeg(data, length, BRUNSLI_READ_ALL, &jpg, nullptr);
  if (status != BRUNSLI_OK) {
    printf("Decoding Brunsli failed with status %d\n", status);
    return result;
  }

  result = new std::string();
  auto write = [](void* data, const uint8_t* buf, size_t count) -> int {
    std::string* result = reinterpret_cast<std::string*>(data);
    result->append(reinterpret_cast<const char*>(buf), count);
    return count;
  };
  JPEGOutput writer(write, result);
  if (!WriteJpeg(jpg, writer)) {
    printf("Serializing JPEG failed\n");
    delete result;
    result = nullptr;
  }
  return result;
}

const char* GetJpegData(std::string* jpeg) {
  if (!jpeg) return 0;
  return jpeg->data();
}

size_t GetJpegLength(std::string* jpeg) {
  if (!jpeg) return 0;
  return jpeg->size();
}

void FreeJpeg(std::string* jpeg) { delete jpeg; }

std::string* JpegToBrunsli(const uint8_t* data, size_t length) {
  std::string* result = nullptr;
  brunsli::JPEGData jpg;
  if (!ReadJpeg(data, length, JPEG_READ_ALL, &jpg)) {
    printf("Parsing JPEG failed\n");
    return result;
  }
  size_t result_length = brunsli::GetMaximumBrunsliEncodedSize(jpg);
  result = new std::string();
  result->resize(result_length);
  bool ok = brunsli::BrunsliEncodeJpeg(
      jpg, reinterpret_cast<uint8_t*>(&result->at(0)), &result_length);
  if (ok) {
    result->resize(result_length);
  } else {
    printf("Encoding Brunsli failed\n");
    delete result;
    result = nullptr;
  }
  return result;
}

const char* GetBrunsliData(std::string* brunsli) {
  if (!brunsli) return 0;
  return brunsli->data();
}

size_t GetBrunsliLength(std::string* brunsli) {
  if (!brunsli) return 0;
  return brunsli->size();
}

void FreeBrunsli(std::string* brunsli) { delete brunsli; }

}  // extern "C"