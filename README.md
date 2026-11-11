# D-E cac

D-E cac (POLYD_Ecac) calculates charge and discharge energy density and efficiency from electric displacement–electric field and polarization–electric field hysteresis loops. It supports Polyk and Radiant TXT input.

Developed by 陈东洋.

## Use D-E cac

- **Windows x64** — Download `POLYD_Ecac-1.1.0-win-x64.exe` from [Release v1.1.0](https://github.com/interestingphd/D-E-cac/releases/tag/v1.1.0). The executable includes its runtime.
- **Python source** — With Python 3.11 and Tk installed, install Matplotlib and run the script:

  ```text
  python -m pip install matplotlib==3.11.0
  python polyd_ecac.py
  ```

Choose your input and output folders, select Polyk or Radiant, and click **Run calculation**. Use **Open CSV** to view existing results. The default plot compares discharged energy density (J/cm³) and efficiency (%) against maximum electric field (MV/m). Export the selected plot as PNG, SVG or PDF.

For Radiant, the first two whitespace-separated numeric columns are electric field (kV/cm) and polarization (µC/cm²). The input folder is scanned recursively, with groups identified by relative folder path.

## Cite

- [10.1002/adma.202402239](https://www.x-mol.com/paperRedirect/1771775839390064640)
- [10.1016/j.cclet.2025.111436](https://www.x-mol.com/paperRedirect/1932955793206231040)

## License

The original D-E cac code and documentation are available under the MIT License. Bundled libraries retain their own licenses. Run `POLYD_Ecac-1.1.0-win-x64.exe --licenses` to view their notices.

Copyright (c) 2026 陈东洋

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
