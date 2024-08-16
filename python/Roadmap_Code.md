project_directory/
│

├── python/

│ ├── final.py # Main script

│ ├── co57_isotope.py # Contains 4 inner classes for 12x12 and 8x8 SiPM arrangement regarding to  <sup>57</sup>Co isotope

│ ├── cs137_isotope.py # Contains 4 inner classes for 12x12 and 8x8 SiPM arrangement regarding to  <sup>137</sup>Cs isotope

│ ├── co60_isotope.py # Contains 4 inner classes for 12x12 and 8x8 SiPM arrangement regarding to  <sup>60</sup>Co isotope

│ ├── scalability.py # Contains classes for future projects

│ ├── config.py # Contain global variables

│ └── final_code_all_isotopes.py # Contains final code in one file



<h2>Detailed Python Code</h2>

<p>
    The Python code [3.8] developed for the accuracy of azimuthal angle determination of the radioactive isotopes <sup>57</sup>Co, <sup>137</sup>Cs, and <sup>60</sup>Co can be found on the GitHub repository named <a href="https://github.com/your-repo-name/gamma-ray-360-detection">"gamma-ray-360-detection"</a>.
</p>

<p>
    The program consists of three main classes, each with four subclasses. The first two classes include the analysis for the 12x12 SiPM array setup, while the last two cover the 8x8 SiPM array setup.
</p>

<h3>1. <sup>57</sup>Co Isotope</h3>
<ul>
    <li><strong>Co57_AnalyticalMethodSolution_12x12_SiPM</strong></li>
    <li><strong>Co57_DataAnalysis_12x12_SiPM</strong></li>
    <li><strong>Co57_AnalyticalMethodSolution_8x8_SiPM</strong></li>
    <li><strong>Co57_DataAnalysis_8x8_SiPM</strong></li>
</ul>

<h3>2. <sup>137</sup>Cs Isotope</h3>
<ul>
    <li><strong>Cs137_AnalyticalMethodSolution_12x12_SiPM</strong></li>
    <li><strong>Cs137_DataAnalysis_12x12_SiPM</strong></li>
    <li><strong>Cs137_AnalyticalMethodSolution_8x8_SiPM</strong></li>
    <li><strong>Cs137_DataAnalysis_8x8_SiPM</strong></li>
</ul>

<h3>3. <sup>60</sup>Co Isotope</h3>
<ul>
    <li><strong>Co60_AnalyticalMethodSolution_12x12_SiPM</strong></li>
    <li><strong>Co60_DataAnalysis_12x12_SiPM</strong></li>
    <li><strong>Co60_AnalyticalMethodSolution_8x8_SiPM</strong></li>
    <li><strong>Co60_DataAnalysis_8x8_SiPM</strong></li>
</ul>

<p>
    The classes "X AnalyticalMethodSolution YxY SiPM" isolate the events of the photopeak for each isotope and calculate the error of the azimuthal angle for each position of the isotope, comparing the actual angle with the angle determined by the center-of-mass vector algorithm. They accept input .txt files from a user-specified path for each angle, for a specific number of events, e.g., (10,000 events).
</p>

<p>
    The classes "X DataAnalysis YxY SiPM" are responsible for displaying the normal Gaussian distribution resulting from the analysis, which is of primary interest to the reader.
</p>

