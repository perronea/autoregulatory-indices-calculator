# Caclulate secondary reactivity indices to assess cerebral autoregulation
reactivityIndices.R takes in raw physiological data (ICP, ABP, MAP, ECG, cerebral blood oxygenation) from Mobius medical medical systems and a Niro 300 near-infrared spectrometer and outputs various reactivity indices (moving correlation coefficients) which can be used to determine if a patient is undergoing impaired cerebral autoregulation.

The purpose of this program is to determine cerebral autoregulation in neurocritical care patients to guide clinical management and predict outcome. Cerebral autoregulation is a vital hemodynamic function of cerebral vasculature designed to maintain adequate blood supply to the brain. Dysautoregulation occurs when a patient undergoes severe TBI, hemorrhagic or ischemic stroke, hydrocephalaus, or other brain injuries that result in edema or hypovolemia. Following such an event the systems that control autoregulation become impaired and the patient may suffer from secondary insults such as delayed cerebral ischemia caused by insufficient cerebral blood flow.

The correlation coefficient is calculated as a simple Pearson correlation coefficient that ranges from +1 to -1. When cerebral autoregulation is effective, a low or negative correlation will be present between the input signal (e.g. CPP) and the output signal (e.g. TCD flow). When cerebral autoregulation is ineffective, the input signal will be transmitted passively to the output signal and this will lead to high correlations.<sup>1</sup>


<sup>1</sup> Zweifel C. Dias C. Smielewski P. Czosnyka M., Continuous time-domain monitoring of cerebral autoregulation in neurocritical care. Medical Engineering and Physics, 2014 vol. 35 (5) pp:638-635
