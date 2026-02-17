# Latex
* New line for each sentence for ease of tracking.
* Use -- to connect two names, e.g., Euler--Lagrange. 
* First sententce in abstract research context (challenge, why we do research). Then, second sentence what we did. 
* Be specific about how much improvement we have, avoid significant etc vague words.
* Math fonts: lower case, regular for scalar. lower case, bold for vector. Upper case bold for matrix. (exception are random variables?)
* Consistency of variables. Avoid overloaded variables.
* When first define variable, make sure define the sizes. When new var show up, define it.
* When ready for coauthor to edit, make sure push all changaes, download a clean version and compiles.
* Avoid empyty section and subsection --  e.g. \section{xxc} directly followed by \subsection{xxx} no words in between.
* Avoid one sentence paragraph.
* Avoid ending a paragraph with equation.
* Avoid the words "given by" unless clear "given by" whom.
* Avoid exp use e^.
* Avoid one subsection section. same logic apply recursively.
* Avoid overfull box.
* Put all newly generated bib in a separate file named agi.bib. 
* In examples, refer to the genberal equation and define parameters.
* Avoid bullet point if not absolutely necessary.
* Check acronyms. Define once and only once. 
* In paragraph eqn frac use xx/xx. Out paragraph eqn opposite.
* Make sure all tables, figures and appendices are refed.
* Use coloneqq instead of :=.




# Figure
* Use niceplot. Math fonts use CMU.
* y axis rotate horizontal.
* Avoid grid etc noise in the figure.
* Sequential data use sequatial color scheme. Diverging data use diverging color scheme (pay attention that zero needs to be white).
* Avoid dashed lines. Use color transparency etc.
* When comparing two lines and showing they are close, we can plot one line wider and more transparent and the other thinner and darker.
* Figure font size similar with the main text.
* Avoid figure subtitle/title in the paper file. Instead if multiple subfigure title needed add to the latex.
* Avoid text in the figure overlap with lines--it is okay to extend the figure a bit to create new space if too crowded.
* For horizontal multiple figures if share y axis and vertical multiple figures the share x axis, only show left y and lower x axis with ticks and labels. All rest axis removed.
* If multuple figure share one color bar only keep one.
* colorbar put it under the figure.
* If the variables are integer value, avoid using 12.4 etc fraction number for ticks.
* Use few ticks--two ticks fully determine the scale (but we dont need to be that extreme).
* If same color labels used across a figure, show the label once is okay. (Ideally in the right top figure).
* Better place the figures, especially avoid many figures after the text that discuss them making the figures shwoing up in other sections.

# Code
* (Recomended) Use a separate config.py for hyperparameters and problem-specific settings. Avoid hardcoding values in main scripts.
* Extract shared utilities (model classes, common functions) into a single file (or class) at the parent directory level (e.g., flowmap.py) to avoid code duplication across examples.
* Use clear script names: generate_data.py, train.py, eval.py (or similar pipeline stages).
* Set random seeds explicitly at the start of each script for reproducibility.
* Keep example directories self-contained: each example should have its own config.py defining problem-specific drift, diffusion, ranges, etc.
