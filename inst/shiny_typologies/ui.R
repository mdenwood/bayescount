library('shiny')
library('shinythemes')

# logifySlider javascript function
JS.logify <-
"
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.round(Math.pow(10, num)*100)/100); }
    })
  }
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('k1', sci = false)
    logifySlider('k2', sci = false)
  }, 5)})
"


# Default min, max and step for log sliders:
logslidvals <- c(-2,2,0.025)

fluidPage(
	# For log scale sliders:
    tags$head(tags$script(HTML(JS.logify))),
    tags$head(tags$script(HTML(JS.onload))),

	theme = shinytheme("cerulean"),
	
	hr(),
#	h3("Select Parameter Values", style="text-align:center; "),
#	hr(),
	
	tabsetPanel(
		tabPanel("Parameters",
			fluidRow(
				column(6, 
			 
					 # div(style = "height:200px;background-color: yellow;", "Topleft"),
					 h4("Choose Study Design Parameters", style="text-align:center; "),
					 hr(),

					 selectInput("type", "Study Type", c("Simple Paired"), selected="Simple Paired", width='100%'),
					 # Also allow Simple Unpaired and more complex designs with replicates, pooling, unequal EDT etc?

					 hr(),

					 numericInput('N', 'Sample Size', value=10, min=5, step=1, width='100%'),
					 numericInput('mu', 'Mean pre-treatment count (EDT * mean EPG)', value=20, min=1, step=0.1, width='100%')
			 
				),
				column(6,
					h4("Choose Efficacy Parameters", style="text-align:center; "),
					hr(),
					numericInput('target', 'Target efficay (%)', min=0, max=100, value=95, step=0.1, width='100%'),
					numericInput('delta', 'Non-inferiority margin (% points)', min=0, max=25, value=5, step=0.1, width='100%'),
					numericInput('pthresh', 'Threshold for significance (p)', min=0, max=0.1, value=0.025, step=0.025, width='100%')
				)
			),
			fluidRow(
				hr(),
				h4("Choose Over-Dispersion Parameters", style="text-align:center; "),
				hr(),

				column(6,
					selectInput("kpreset", "Preset k values", c("Sheep", "Cattle", "Calves", "Equine", "Custom"), selected="Sheep", width='100%'),
	 
					conditionalPanel(
				 		condition = "input.kpreset == 'Custom'",
						sliderInput('k1', 'Pre-treatment k', min=logslidvals[1], max=logslidvals[2], value=0, step=logslidvals[3], width='100%'),
						sliderInput('k2', 'Post-treatment k', min=logslidvals[1], max=logslidvals[2], value=-0.30, step=logslidvals[3], width='100%'),
						sliderInput('cor', 'Correlation', min=0.01, max=0.99, value=0.2, step=0.01, width='100%')
					)
				),
				column(6,
					htmlOutput('kvalues')
				)
			),

			hr(),
			actionButton("simulate", "Click to run the simulation with these parameters", width="100%"),
			hr(),
			htmlOutput('footer'),
			hr()
		),
		
		tabPanel("Probabilities",
			fluidRow(
				column(6,
					h4("Choose Statistical Method", style="text-align:center; "),
					selectInput("statmethod", "Statistical Method", c("BNB", "WAAVP", "Gamma", "Dobson", "Asymptotic"), selected="BNB", width='100%')
				),
				column(6,
					h4("Choose Statistical Framework", style="text-align:center; "),
					selectInput("statframe", "Statistical Framework", c("Denwood", "Kaplan", "Coles"), selected="Denwood", width='100%')
				)
			),
			fluidRow(
				hr(),
				htmlOutput('parameters'),
				hr(),
				h4("Probability of a non-inconclusive result"),
				plotOutput("noninc_plot"),
				hr(),
				h4("Probability of each classification type"),
				plotOutput("typologies_plot")
			)
		),
		
		tabPanel("Typologies",
			fluidRow(
				h4("Description of typologies", style="text-align:center; "),
				column(6,
					h4("Typologies"),
					imageOutput("typologies")
				),
				column(6,
					h4("Classifications"),
					htmlOutput('classifications')
				)
			)
		)
	
	)
)

