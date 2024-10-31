# Changelog

## 0.0.2 (2024-10-31)
1. Python Changes
    - Removed JSON support completely
    - Parameter dataclass can not be addressed as dictionaries
    - Internally, `HMSPredictor.sigma` is still groove separation
2. Configuration Format Changes
    - Removed `SysParam.hmsVersion`
    - Changed `SysParam.sigma` to grating groove density instead of groove separation