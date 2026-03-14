#' Initialise hCoCena Object
#' 
#' Creates an object called 'hcobject' in the global envirentment that will be used by hCoCena throughout the analysis.
#' @param envo Environment where `hcobject` should be created. Defaults to `.GlobalEnv`.
#' @export

init_object <- function(envo = .GlobalEnv){
	.hc_legacy_warning("init_object")
	hcobject <- as_hcobject(hc_init())
	.hc_assign_legacy_hcobject(hcobject)
	if(!base::identical(envo, .hc_legacy_state_env())){
		base::assign("hcobject", hcobject, envir = envo)
	}
}
