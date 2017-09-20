# na wejscie obiekt (Expression Set), na wyjscie ten sam obiekt z podmienionym cdf-em

custom_cdf = function(obiekt){

  # nazwa cdf dla danego obiektu
  cdf_name = cleancdfname(cdfName(obiekt), addcdf = FALSE)
  #dlugosc nazwy
  dl_cdf = nchar(cdf_name)
  #system operacyjny
  sys_op = Sys.info()['sysname']

  # sprawdzenie czy nazwa zaczyna się od hgu lub hth (w bazie BrainArray biblioteki te pisane wielkimi literami)
  if (substr(cdf_name, 1, 3) == 'hgu' || substr(cdf_name, 1, 3) == 'hth'){

    cdf_name = toupper(cdf_name)
    # w razie, gdyby w bibliotece wystąpiło słowo plus (ze względu na nazwy bibliotek w BrainArray)
    cdf_name = sub("PLUS" , "Plus", cdf_name, ignore.case = FALSE, perl = FALSE)
  }

  # w Brain Array nie ma informacji o wersji poza tymi w if-ie
  if (substr(cdf_name, 1, 6) != 'HGU95A' && substr(cdf_name, 1, 8) != 'u133aaofa'){
    if (substr(cdf_name, dl_cdf-1, dl_cdf-1) == 'v'){
      cdf_name = substr(cdf_name, 1, dl_cdf-2)
    }
  }

  # nazwa pliku do pobrania
  nazwa_pobieranego_pliku = paste(tolower(cdf_name), "hsentrezgcdf", sep = "")

  # jeśli nie ma jeszcze takiej paczki -> pobranie, instalacja i wczytanie;
  # jeśli paczka już jest zainstalowana -> wczytanie
  if(!(nazwa_pobieranego_pliku %in% row.names(installed.packages()))){
    message('Please wait a bit while annotation package is being downloaded, installed and loaded.')

    # if(sys_op=='Linux'){
    download.file(paste('http://mbni.org/customcdf/20.0.0/entrezg.download/', nazwa_pobieranego_pliku, '_20.0.0.tar.gz', sep = ''), paste(getwd(), '/', nazwa_pobieranego_pliku, '_20.0.0.tar.gz', sep = ''))
    install.packages(paste(getwd(), '/', nazwa_pobieranego_pliku, '_20.0.0.tar.gz', sep = ''), type='source', repos = NULL)
    # }else if(sys_op=='Windows'){
    #   download.file(paste('http://mbni.org/customcdf/20.0.0/entrezg.download/', nazwa_pobieranego_pliku, '_20.0.0.zip', sep = ''), paste(getwd(), '/', nazwa_pobieranego_pliku, '_20.0.0.tar.gz', sep = ''))
    #   install.packages(paste(getwd(), '/', nazwa_pobieranego_pliku, '_20.0.0.zip', sep = ''), repos = NULL)
    # }
    library(nazwa_pobieranego_pliku, character.only = TRUE)

  }else{
    message('Please wait a bit while annotation package is being loaded.')
    library(nazwa_pobieranego_pliku, character.only = TRUE)
  }

  obiekt@cdfName = paste(cdf_name, '_Hs_ENTREZG', sep = '')

  return(obiekt)
}
