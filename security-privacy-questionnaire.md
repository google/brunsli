Document outline copied from https://www.w3.org/TR/security-privacy-questionnaire/

### 2. Questions to Consider

##### 2.1. What information might this feature expose to Web sites or other parties, and for what purposes is that exposure necessary?
Content-Encoding does not affect the information exposure.

##### 2.2. Is this specification exposing the minimum amount of information necessary to power the feature?
Yes.

##### 2.3. How does this specification deal with personal information or personally-identifiable information or information derived thereof?
Content-Encoding does not affect the way personal information, PII and derivatives are processed.

##### 2.4. How does this specification deal with sensitive information?
Content-Encoding does not affect the way sensitive information is processed.

##### 2.5. Does this specification introduce new state for an origin that persists across browsing sessions?
No.

##### 2.6. What information from the underlying platform, e.g. configuration data, is exposed by this specification to an origin?
None.

##### 2.7. Does this specification allow an origin access to sensors on a user’s device
No.

##### 2.8. What data does this specification expose to an origin? Please also document what data is identical to data exposed by other features, in the same or different contexts.
Only extra "jxl" item in "Accept-Encodings" HTTP request header.

##### 2.9. Does this specification enable new script execution/loading mechanisms?
No.

##### 2.10. Does this specification allow an origin to access other devices?
No.

##### 2.11. Does this specification allow an origin some measure of control over a user agent’s native UI?
No.

##### 2.12. What temporary identifiers might this this specification create or expose to the web?
None.

##### 2.13. How does this specification distinguish between behavior in first-party and third-party contexts?
Not distinguished.

##### 2.14. How does this specification work in the context of a user agent’s Private Browsing or "incognito" mode?
Same way in regular context.

##### 2.15. Does this specification have a "Security Considerations" and "Privacy Considerations" section?
Security consideration: parses data.
Mitigation: has been fuzzed for many months.

Security consideration: potential "zip-bomb".
Mitigation: streaming processing + limiting output.

No privacy considerations.

##### 2.16. Does this specification allow downgrading default security characteristics?
No.

##### 2.17. What should this questionnaire have asked?
Nothing at the moment.

### 3. Threat Models

##### 3.1. Passive Network Attackers
Compression ratio for most non-nonsenical JPEGs is quite stable over the range of images, and depends only on image quality and encoding parameters.
Not considered as sensible information.

##### 3.2. Active Network Attackers
Attacker might try to implement "zip-bomb" attack to consume victims browser memory.

##### 3.3. Same-Origin Policy Violations
N/A

##### 3.4. Third-Party Tracking
N/A

##### 3.5. Legitimate Misuse
N/A
