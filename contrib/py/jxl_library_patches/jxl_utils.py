# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.
"""Utils for patching JPEG-XL support into Python image libraries.

"""
import os
import shutil
import subprocess
import tempfile


# Header for JPEG-XL's recompressed JPEG mode.
JPEGXL_RECOMPRESSED_JPEG_HEADER = b'\x0a\x04\x42\xd2\xd5\x4e'


# Default binary to use to convert JXL recompressed JPEG to JPEG.
# This is expected to have signature:
#   {jpeg_from_jxl} <infile:jxl> <outfile:jpeg>
#
# _BIN_JPEG_FROM_JXL = 'jpeg_from_jxl'
_BIN_JPEG_FROM_JXL = 'dbrunsli'


def get_bin_jpeg_from_jxl():
  """Gets the JPEG-XL-to-JPEG converter command."""
  return os.getenv('BIN_JPEG_FROM_JXL', _BIN_JPEG_FROM_JXL)


def can_call_converter():
  try:
    subprocess.call(
        [get_bin_jpeg_from_jxl()],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  except OSError:  # Could not call converter.
    return False
  except subprocess.CalledProcessError:
    return True
  return True


def is_jpegxl_recompressed_jpeg_file(filename):
  """Returns True iff the given filename is a genuine JPEG-XL file."""
  try:
    with open(filename, 'rb') as h:
      header = h.read(len(JPEGXL_RECOMPRESSED_JPEG_HEADER))
      # Cf. https://arxiv.org/pdf/1908.03565.pdf, section 9.1,
      # on recompressed-JPEG header.
      return header == JPEGXL_RECOMPRESSED_JPEG_HEADER
  except:  # pylint:disable=bare-except
    # If anything failed, this means that we cannot establish that the file
    # has the expected header, so we return False.
    return False


def call_with_jxl_filename_arg1_replaced_by_temp_jpeg_file(func,
                                                           *args,
                                                           **kwargs):
  """Calls 'func' with the 1st arg from `args` changed to point to a JPEG.

  The JPEG file is a temporary file generated on-the-fly by converting the
  file given by the original 1st argument.

  Args:
    func: The function to call.
    *args: The arguments to use for the function call, with args[0] being the
       JPEG-XL image file name which gets replaced with a JPEG filename for the
       call.
    **kwargs: Passthrough keyword args for `func`.

  Returns:
    The result of calling func(jpeg_filename, *args[1:], **kwargs)

  Raises:
    RuntimeError: if behind-the-scenes image conversion failed.
  """
  # Give the converter some nice workspace where it also can place its own
  # temporary files if it wanted to, without having to worry about security
  # issues.
  tempdir = tempfile.mkdtemp()
  try:
    temp_jpeg_filename = os.path.join(tempdir, 'temp_converted_jxl.jpeg')
    if 0 != subprocess.call(
        [get_bin_jpeg_from_jxl(), args[0], temp_jpeg_filename]):
      raise RuntimeError(
          'Ad-hoc conversion of JPEG-XL to JPEG failed on file: %r' % args[0])
    return func(*((temp_jpeg_filename,) + args[1:]),
                **kwargs)
  finally:
    # Clean up the securely-created temporary directory and everything inside.
    shutil.rmtree(tempdir)
