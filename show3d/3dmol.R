#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= r3dmol -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

mol3d<- function(pdb_view, ligid, check_style){
#browser()
  library(r3dmol)
  # source("bib/common-v2.R")
  # source("bib/fun-base-v46.R")
  # source("bib/fun-mass-v46.R")
  # source("main/fun-ligs-v1.R")
  # source("main/main-ligs.R")
  # 

if(check_style == "spheres") {
  r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      antialias = TRUE,
      lowerZoomLimit = 50,
      upperZoomLimit = 400,
    ),

  )%>%
    
    # Add model to scene
   m_add_model(data = paste0("pdb-mol/", pdb_view), format = "pdb") %>%
   # m_add_model(data = paste0("pdb-mol/", "6XQU.pdb"), format = "pdb") %>%

    # Zoom to encompass the whole scene
    m_zoom_to() %>%
    # m_set_style(style = m_style_stick(radius = 0.1)) %>% #<<<<<<

    m_set_style(style = m_style_sphere(radius = 1.5))%>%
    m_add_res_labels(style = m_style_label(
      
      fontSize = 12,
      fontColor = "Blue",
      fontOpacity = 0.9,
      showBackground = F
    )) %>%
    m_set_style(
      sel = m_sel(resi = ligid), #resn = "ALD" elem="HETATM"
      style = m_style_stick(
        radius = 0.3,
        colorScheme = "magenta"
      )) %>%
    m_set_style(
      sel = m_sel(resn = "HOH"), #resn = "ALD" elem="HETATM"
      style = m_style_sphere(radius = 0.8))%>%
    m_add_surface(style = m_style_surface(opacity = 0.4, colorScheme = "ssPyMOL"), focus = 1 ) %>%
    m_spin(speed = 0)
  }
  else {
    r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        antialias = TRUE,
        lowerZoomLimit = 50,
        upperZoomLimit = 400,
      ),
      
    )%>%
      # Add model to scene
      m_add_model(data = paste0("pdb-mol/", pdb_view), format = "pdb") %>%
      # Zoom to encompass the whole scene
      m_zoom_to() %>%
    m_set_style(style = m_style_stick(radius = 0.1))%>%
    m_add_res_labels(style = m_style_label(
      fontSize = 10,
      fontColor = "Blue",
      fontOpacity = 0.7,
      showBackground = F
    )) %>%
    m_set_style(
      sel = m_sel(resn = lig.table$lig1n),
        # sel = m_sel(resi = ligid), #resn = "ALD" elem="HETATM"
        style = m_style_stick(
          radius = 0.3,
          colorScheme = "magenta"
     )) %>%
    m_set_style(
        sel = m_sel(resn = "HOH"), #resn = "ALD" elem="HETATM"
        style = m_style_sphere(radius = 0.8))%>%
    m_add_surface(style = m_style_surface(opacity = 0.4, colorScheme = "ssPyMOL"), focus = 1 ) %>%
    m_spin(speed = 0)
    }


}

