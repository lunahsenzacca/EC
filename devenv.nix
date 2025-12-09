{ pkgs, lib, config, inputs, ... }:

{
  packages = with pkgs; [
    gcc-unwrapped
    libz
    jdk21_headless
    glib
  ];
  languages.python.enable = true;
  languages.python.version = "3.13.9";
  languages.python.venv.enable = true;
  languages.python.poetry.enable = true;
  languages.python.poetry.install.enable = true;
}
